//! A logger that can hold a thread's lines back and write them later, so the
//! samples of a cohort can run side by side and still log one after another.
//!
//! `--sample all` used to annotate one sample after the next, and every line a
//! sample logged came before every line of the sample after it. Run side by
//! side they would interleave, and most lines do not say which sample they
//! belong to. So each sample's thread holds its lines, and the run writes them
//! in sample order once the samples before have been written.
//!
//! Holding lines needs the process logger to be the one [`install`] sets up.
//! When it is not, a cohort runs its samples one at a time as it always did.

use log::{Level, LevelFilter, Log, Metadata, Record, SetLoggerError};
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, Ordering};

/// Set once [`install`] has made the holding logger the process logger.
static INSTALLED: AtomicBool = AtomicBool::new(false);

thread_local! {
    /// What this thread is holding, innermost hold last. Holds nest: loading
    /// an input the cohort shares keeps that load's lines apart from the sample
    /// that happened to ask first, so every sample can be given them.
    static HELD: RefCell<Vec<Vec<Held>>> = const { RefCell::new(Vec::new()) };
}

/// Something written while a thread was holding: a log line, or what an
/// outside program run meanwhile wrote, which belongs among the same lines.
#[derive(Debug, Clone)]
pub(crate) enum Held {
    Line(LogLine),
    Output { stdout: Vec<u8>, stderr: Vec<u8> },
}

impl Held {
    /// Write it now, or hold it again if this thread is holding.
    pub(crate) fn write(&self) {
        match self {
            Held::Line(line) => line.write(),
            Held::Output { stdout, stderr } => {
                if !hold(|| self.clone()) {
                    use std::io::Write;
                    let _ = std::io::stdout().write_all(stdout);
                    let _ = std::io::stderr().write_all(stderr);
                }
            }
        }
    }
}

/// Put what `make` gives into the innermost hold of this thread, if it is
/// holding, and say whether it was.
fn hold(make: impl FnOnce() -> Held) -> bool {
    HELD.try_with(|held| {
        let holding = !held.borrow().is_empty();
        if holding {
            // Made before the stack is borrowed, in case making it logs.
            let item = make();
            if let Some(innermost) = held.borrow_mut().last_mut() {
                innermost.push(item);
            }
        }
        holding
    })
    .unwrap_or(false)
}

/// One log line, kept until it can be written.
#[derive(Debug, Clone)]
pub(crate) struct LogLine {
    level: Level,
    target: String,
    module_path: Option<String>,
    file: Option<String>,
    line: Option<u32>,
    message: String,
}

impl LogLine {
    fn of(record: &Record<'_>) -> Self {
        LogLine {
            level: record.level(),
            target: record.target().to_string(),
            module_path: record.module_path().map(str::to_string),
            file: record.file().map(str::to_string),
            line: record.line(),
            message: record.args().to_string(),
        }
    }

    /// Log this line again, as it was first logged: held if this thread is
    /// holding lines, written otherwise.
    pub(crate) fn write(&self) {
        log::logger().log(
            &Record::builder()
                .args(format_args!("{}", self.message))
                .level(self.level)
                .target(&self.target)
                .module_path(self.module_path.as_deref())
                .file(self.file.as_deref())
                .line(self.line)
                .build(),
        );
    }
}

struct HoldingLogger<L> {
    inner: L,
}

impl<L: Log> Log for HoldingLogger<L> {
    fn enabled(&self, metadata: &Metadata<'_>) -> bool {
        self.inner.enabled(metadata)
    }

    fn log(&self, record: &Record<'_>) {
        // A line the logger would drop is not worth holding.
        if !self.inner.enabled(record.metadata()) {
            return;
        }
        if !hold(|| Held::Line(LogLine::of(record))) {
            self.inner.log(record);
        }
    }

    fn flush(&self) {
        self.inner.flush();
    }
}

/// Make `inner` the process logger, wrapped so that lines can be held back.
/// Called once, by the program, in place of the logger's own `init`.
pub fn install<L: Log + 'static>(inner: L, max_level: LevelFilter) -> Result<(), SetLoggerError> {
    log::set_boxed_logger(Box::new(HoldingLogger { inner }))?;
    log::set_max_level(max_level);
    INSTALLED.store(true, Ordering::Release);
    Ok(())
}

/// Whether lines can be held back, which takes the logger [`install`] sets up.
pub(crate) fn can_hold() -> bool {
    INSTALLED.load(Ordering::Acquire)
}

/// Run `work` holding back every line it logs on this thread, and hand those
/// lines back with its result, in the order they were logged.
pub(crate) fn holding<R>(work: impl FnOnce() -> R) -> (R, Vec<Held>) {
    /// Ends the hold even if `work` panics, so the thread does not go on
    /// holding lines nobody will write.
    struct Hold;
    impl Drop for Hold {
        fn drop(&mut self) {
            let _ = HELD.try_with(|held| held.borrow_mut().pop());
        }
    }

    HELD.with(|held| held.borrow_mut().push(Vec::new()));
    let hold = Hold;
    let result = work();
    let lines = HELD.with(|held| {
        held.borrow_mut()
            .last_mut()
            .map(std::mem::take)
            .unwrap_or_default()
    });
    drop(hold);
    (result, lines)
}

/// How many lines the innermost hold of this thread has kept so far, which is
/// where anything reported now falls among them.
pub(crate) fn held_so_far() -> usize {
    HELD.with(|held| held.borrow().last().map_or(0, Vec::len))
}

/// Run an outside program. While this thread holds its lines, what the
/// program writes is held with them, so it stays among the lines of the
/// sample that ran it; otherwise the program writes where this process does.
pub(crate) fn run_program(
    command: &mut std::process::Command,
) -> std::io::Result<std::process::ExitStatus> {
    let holding = HELD
        .try_with(|held| !held.borrow().is_empty())
        .unwrap_or(false);
    if !holding {
        return command.status();
    }
    let output = command.output()?;
    hold(|| Held::Output {
        stdout: output.stdout,
        stderr: output.stderr,
    });
    Ok(output.status)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::sync::{Mutex, OnceLock};

    /// What the logger behind the holding one received, in order, with the
    /// thread that wrote it: the tests of the crate run at once, and a test
    /// tells its own lines apart by its thread.
    fn written() -> &'static Mutex<Vec<(std::thread::ThreadId, String)>> {
        static WRITTEN: OnceLock<Mutex<Vec<(std::thread::ThreadId, String)>>> = OnceLock::new();
        WRITTEN.get_or_init(|| Mutex::new(Vec::new()))
    }

    /// Every line this thread has written so far.
    pub(crate) fn written_here() -> Vec<String> {
        let here = std::thread::current().id();
        written()
            .lock()
            .unwrap()
            .iter()
            .filter(|(thread, _)| *thread == here)
            .map(|(_, line)| line.clone())
            .collect()
    }

    struct Recorder;

    impl Log for Recorder {
        fn enabled(&self, metadata: &Metadata<'_>) -> bool {
            metadata.level() <= Level::Info
        }
        fn log(&self, record: &Record<'_>) {
            if self.enabled(record.metadata()) {
                written().lock().unwrap().push((
                    std::thread::current().id(),
                    format!("{} {}", record.target(), record.args()),
                ));
            }
        }
        fn flush(&self) {}
    }

    /// The tests of this crate share one process logger, so whichever needs
    /// it first installs it for all of them.
    pub(crate) fn install_for_tests() {
        static ONCE: OnceLock<()> = OnceLock::new();
        ONCE.get_or_init(|| {
            install(Recorder, LevelFilter::Info).expect("no other logger in the unit tests");
        });
    }

    fn written_with(marker: &str) -> Vec<String> {
        written()
            .lock()
            .unwrap()
            .iter()
            .filter(|(_, line)| line.contains(marker))
            .map(|(_, line)| line.clone())
            .collect()
    }

    #[test]
    fn held_lines_are_written_when_asked_and_in_order() {
        install_for_tests();
        let ((), lines) = holding(|| {
            log::info!(target: "held-test", "first");
            log::debug!(target: "held-test", "below the level, never kept");
            log::info!(target: "held-test", "second");
        });
        assert!(written_with("held-test").is_empty(), "nothing written yet");
        assert_eq!(lines.len(), 2);
        for line in &lines {
            line.write();
        }
        assert_eq!(
            written_with("held-test"),
            vec!["held-test first", "held-test second"]
        );
    }

    /// A line logged while an inner hold is open goes to that hold, and what
    /// the inner hold gave back can be logged into the outer one.
    #[test]
    fn holds_nest() {
        install_for_tests();
        let ((), outer) = holding(|| {
            log::info!(target: "nest-test", "a");
            let ((), inner) = holding(|| log::info!(target: "nest-test", "b"));
            assert_eq!(held_so_far(), 1);
            for line in &inner {
                line.write();
            }
            assert_eq!(held_so_far(), 2);
            log::info!(target: "nest-test", "c");
        });
        let messages: Vec<&str> = outer
            .iter()
            .map(|held| match held {
                Held::Line(line) => line.message.as_str(),
                Held::Output { .. } => "output",
            })
            .collect();
        assert_eq!(messages, ["a", "b", "c"]);
        assert!(written_with("nest-test").is_empty());
    }

    #[test]
    fn a_thread_that_holds_nothing_writes_at_once() {
        install_for_tests();
        std::thread::spawn(|| log::info!(target: "direct-test", "now"))
            .join()
            .unwrap();
        assert_eq!(written_with("direct-test"), vec!["direct-test now"]);
    }
}
