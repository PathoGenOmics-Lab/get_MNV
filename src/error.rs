//! Structured error types with categorised exit codes and optional source chains.

use serde::Serialize;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs::File;
use std::io::Write;

/// Error category, each with a unique mnemonic code and process exit code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorCode {
    Generic,
    Config,
    Validation,
    Io,
    Csv,
    Htslib,
    Utf8,
    ParseInt,
    ParseFloat,
    Json,
}

impl ErrorCode {
    pub fn as_str(self) -> &'static str {
        match self {
            ErrorCode::Generic => "E000",
            ErrorCode::Config => "E001",
            ErrorCode::Validation => "E002",
            ErrorCode::Io => "E100",
            ErrorCode::Csv => "E101",
            ErrorCode::Htslib => "E102",
            ErrorCode::Utf8 => "E103",
            ErrorCode::ParseInt => "E104",
            ErrorCode::ParseFloat => "E105",
            ErrorCode::Json => "E106",
        }
    }

    pub fn exit_code(self) -> i32 {
        match self {
            ErrorCode::Generic => 1,
            ErrorCode::Config => 2,
            ErrorCode::Validation => 3,
            ErrorCode::Io => 10,
            ErrorCode::Csv => 11,
            ErrorCode::Htslib => 12,
            ErrorCode::Utf8 => 13,
            ErrorCode::ParseInt => 14,
            ErrorCode::ParseFloat => 15,
            ErrorCode::Json => 16,
        }
    }
}

impl Display for ErrorCode {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.write_str(self.as_str())
    }
}

/// Application error with a category code, human message, and optional cause.
#[derive(Debug)]
pub struct AppError {
    pub code: ErrorCode,
    pub message: String,
    source: Option<Box<dyn std::error::Error + Send + Sync>>,
}

impl Display for AppError {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "[{}] {}", self.code, self.message)?;
        if let Some(ref cause) = self.source {
            write!(f, ": {}", cause)?;
        }
        Ok(())
    }
}

impl std::error::Error for AppError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_ref()
            .map(|e| e.as_ref() as &(dyn std::error::Error + 'static))
    }
}

impl AppError {
    /// Create an error with a specific code and message.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::error::{AppError, ErrorCode};
    /// let err = AppError::new(ErrorCode::Validation, "invalid allele");
    /// assert_eq!(err.code, ErrorCode::Validation);
    /// ```
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            source: None,
        }
    }

    /// Attach an underlying cause to this error.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::error::AppError;
    /// let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "gone");
    /// let err = AppError::io("open failed").with_source(io_err);
    /// assert!(err.to_string().contains("gone"));
    /// ```
    pub fn with_source(
        mut self,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        self.source = Some(Box::new(source));
        self
    }

    /// Add context to an existing error, wrapping it as the source.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::error::AppError;
    /// let err = AppError::validation("bad allele");
    /// let wrapped = err.context("while parsing VCF line 42");
    /// assert!(wrapped.to_string().contains("line 42"));
    /// ```
    pub fn context(self, ctx: impl Into<String>) -> Self {
        let code = self.code;
        let ctx_msg = ctx.into();
        AppError {
            code,
            message: ctx_msg,
            source: Some(Box::new(self)),
        }
    }

    pub fn msg(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Generic, message)
    }

    pub fn config(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Config, message)
    }

    pub fn validation(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Validation, message)
    }

    pub fn io(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Io, message)
    }

    pub fn json(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Json, message)
    }

    /// Return the full causal chain as a multi-line string (for debugging).
    pub fn chain(&self) -> String {
        let mut lines = vec![self.to_string()];
        let mut current: &dyn std::error::Error = self;
        while let Some(cause) = current.source() {
            lines.push(format!("  caused by: {}", cause));
            current = cause;
        }
        lines.join("\n")
    }
}

// --- From conversions (all preserve the original error as source) ---

impl From<&str> for AppError {
    fn from(value: &str) -> Self {
        AppError::msg(value)
    }
}

impl From<String> for AppError {
    fn from(value: String) -> Self {
        AppError::msg(value)
    }
}

impl From<std::io::Error> for AppError {
    fn from(value: std::io::Error) -> Self {
        let msg = format!("I/O error: {}", value);
        AppError::new(ErrorCode::Io, msg).with_source(value)
    }
}

impl From<csv::Error> for AppError {
    fn from(value: csv::Error) -> Self {
        let msg = format!("CSV error: {}", value);
        AppError::new(ErrorCode::Csv, msg).with_source(value)
    }
}

impl From<rust_htslib::errors::Error> for AppError {
    fn from(value: rust_htslib::errors::Error) -> Self {
        let msg = format!("HTSlib error: {}", value);
        AppError::new(ErrorCode::Htslib, msg).with_source(value)
    }
}

impl From<std::str::Utf8Error> for AppError {
    fn from(value: std::str::Utf8Error) -> Self {
        let msg = format!("UTF-8 error: {}", value);
        AppError::new(ErrorCode::Utf8, msg).with_source(value)
    }
}

impl From<std::string::FromUtf8Error> for AppError {
    fn from(value: std::string::FromUtf8Error) -> Self {
        let msg = format!("UTF-8 conversion error: {}", value);
        AppError::new(ErrorCode::Utf8, msg).with_source(value)
    }
}

impl From<std::num::ParseIntError> for AppError {
    fn from(value: std::num::ParseIntError) -> Self {
        let msg = format!("Integer parsing error: {}", value);
        AppError::new(ErrorCode::ParseInt, msg).with_source(value)
    }
}

impl From<std::num::ParseFloatError> for AppError {
    fn from(value: std::num::ParseFloatError) -> Self {
        let msg = format!("Float parsing error: {}", value);
        AppError::new(ErrorCode::ParseFloat, msg).with_source(value)
    }
}

impl From<serde_json::Error> for AppError {
    fn from(value: serde_json::Error) -> Self {
        let msg = format!("JSON error: {}", value);
        AppError::new(ErrorCode::Json, msg).with_source(value)
    }
}

impl From<std::time::SystemTimeError> for AppError {
    fn from(value: std::time::SystemTimeError) -> Self {
        let msg = format!("System time error: {}", value);
        AppError::new(ErrorCode::Generic, msg).with_source(value)
    }
}

pub type AppResult<T> = Result<T, AppError>;

/// Extension trait to add file-path context to `io::Result`.
pub trait IoResultExt<T> {
    /// Wrap an I/O error with the path that caused it.
    fn with_path(self, path: &str) -> AppResult<T>;
}

impl<T> IoResultExt<T> for std::io::Result<T> {
    fn with_path(self, path: &str) -> AppResult<T> {
        self.map_err(|e| {
            let msg = format!("I/O error on '{}': {}", path, e);
            AppError::new(ErrorCode::Io, msg).with_source(e)
        })
    }
}

/// Extension trait to add context to any `AppResult`.
pub trait ResultExt<T> {
    /// Wrap the error with additional context while preserving the original cause.
    fn context(self, ctx: impl Into<String>) -> AppResult<T>;
}

impl<T> ResultExt<T> for AppResult<T> {
    fn context(self, ctx: impl Into<String>) -> AppResult<T> {
        self.map_err(|e| e.context(ctx))
    }
}

// --- JSON serialisation for structured error output ---

#[derive(Debug, Serialize)]
struct ErrorJson<'a> {
    schema_version: &'static str,
    code: &'a str,
    exit_code: i32,
    message: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    cause: Option<String>,
}

pub fn error_to_json(error: &AppError) -> String {
    let cause = error
        .source
        .as_ref()
        .map(|s| s.to_string());
    let payload = ErrorJson {
        schema_version: "1.0.0",
        code: error.code.as_str(),
        exit_code: error.code.exit_code(),
        message: &error.message,
        cause,
    };
    serde_json::to_string_pretty(&payload)
        .unwrap_or_else(|_| "{\"schema_version\":\"1.0.0\",\"code\":\"E000\",\"exit_code\":1,\"message\":\"Failed to serialize error\"}".to_string())
        + "\n"
}

pub fn write_error_json(path: &str, error: &AppError) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    let cause = error
        .source
        .as_ref()
        .map(|s| s.to_string());
    let payload = ErrorJson {
        schema_version: "1.0.0",
        code: error.code.as_str(),
        exit_code: error.code.exit_code(),
        message: &error.message,
        cause,
    };
    serde_json::to_writer_pretty(&mut file, &payload)?;
    file.write_all(b"\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_code_formatting() {
        assert_eq!(ErrorCode::Config.to_string(), "E001");
        assert_eq!(ErrorCode::Htslib.to_string(), "E102");
        assert_eq!(ErrorCode::Validation.exit_code(), 3);
    }

    #[test]
    fn test_error_message_includes_code() {
        let err = AppError::validation("invalid interval");
        assert!(err.to_string().starts_with("[E002]"));
    }

    #[test]
    fn test_new_error_codes() {
        assert_eq!(ErrorCode::ParseFloat.as_str(), "E105");
        assert_eq!(ErrorCode::Json.as_str(), "E106");
        assert_eq!(ErrorCode::ParseFloat.exit_code(), 15);
        assert_eq!(ErrorCode::Json.exit_code(), 16);
    }

    #[test]
    fn test_with_source_preserves_chain() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file gone");
        let app_err = AppError::io("Failed to open config").with_source(io_err);
        assert!(app_err.to_string().contains("Failed to open config"));
        assert!(app_err.to_string().contains("file gone"));
        assert!(std::error::Error::source(&app_err).is_some());
    }

    #[test]
    fn test_context_wraps_original() {
        let inner = AppError::validation("bad allele");
        let outer = inner.context("while parsing VCF line 42");
        assert!(outer.to_string().contains("while parsing VCF line 42"));
        assert!(outer.chain().contains("bad allele"));
        assert_eq!(outer.code, ErrorCode::Validation);
    }

    #[test]
    fn test_chain_format() {
        let io_err = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "access denied");
        let app_err = AppError::io("cannot write output").with_source(io_err);
        let wrapped = app_err.context("pipeline failed");
        let chain = wrapped.chain();
        assert!(chain.contains("pipeline failed"));
        assert!(chain.contains("caused by:"));
        assert!(chain.contains("cannot write output"));
        assert!(chain.contains("access denied"));
    }

    #[test]
    fn test_io_result_ext_with_path() {
        let result: std::io::Result<()> = Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "no such file",
        ));
        let app_result = result.with_path("/tmp/missing.vcf");
        let err = app_result.unwrap_err();
        assert!(err.message.contains("/tmp/missing.vcf"));
        assert_eq!(err.code, ErrorCode::Io);
    }

    #[test]
    fn test_result_ext_context() {
        let result: AppResult<()> = Err(AppError::validation("bad input"));
        let contexted = result.context("processing sample X");
        let err = contexted.unwrap_err();
        assert!(err.to_string().contains("processing sample X"));
        assert!(err.chain().contains("bad input"));
    }

    #[test]
    fn test_from_parse_float_error() {
        let err: AppError = "not_a_float".parse::<f64>().unwrap_err().into();
        assert_eq!(err.code, ErrorCode::ParseFloat);
        assert!(err.to_string().contains("Float parsing error"));
    }

    #[test]
    fn test_from_serde_json_error() {
        let bad_json = serde_json::from_str::<serde_json::Value>("{bad}");
        let err: AppError = bad_json.unwrap_err().into();
        assert_eq!(err.code, ErrorCode::Json);
        assert!(err.to_string().contains("JSON error"));
    }

    #[test]
    fn test_error_json_includes_cause() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file gone");
        let app_err = AppError::io("Failed to open").with_source(io_err);
        let json = error_to_json(&app_err);
        assert!(json.contains("\"cause\""));
        assert!(json.contains("file gone"));
    }

    #[test]
    fn test_error_json_omits_cause_when_none() {
        let app_err = AppError::validation("bad input");
        let json = error_to_json(&app_err);
        assert!(!json.contains("cause"));
    }

    #[test]
    fn test_all_error_codes_unique() {
        let codes = [
            ErrorCode::Generic,
            ErrorCode::Config,
            ErrorCode::Validation,
            ErrorCode::Io,
            ErrorCode::Csv,
            ErrorCode::Htslib,
            ErrorCode::Utf8,
            ErrorCode::ParseInt,
            ErrorCode::ParseFloat,
            ErrorCode::Json,
        ];
        let strs: Vec<_> = codes.iter().map(|c| c.as_str()).collect();
        let exits: Vec<_> = codes.iter().map(|c| c.exit_code()).collect();
        let mut unique_strs = strs.clone();
        unique_strs.sort();
        unique_strs.dedup();
        assert_eq!(strs.len(), unique_strs.len(), "duplicate error code strings");
        let mut unique_exits = exits.clone();
        unique_exits.sort();
        unique_exits.dedup();
        assert_eq!(exits.len(), unique_exits.len(), "duplicate exit codes");
    }
}
