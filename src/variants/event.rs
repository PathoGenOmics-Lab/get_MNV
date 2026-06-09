//! Allele-level event decomposition shared by VCF/iVar parsing, annotation and
//! output.  The goal is to keep one canonical interpretation of `REF`/`ALT`
//! while preserving the original VCF-compatible alleles.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleEventClass {
    Reference,
    Snp,
    Mnv,
    Insertion,
    Deletion,
    Delins,
    ComplexIndel,
    Symbolic,
}

impl AlleleEventClass {
    pub fn as_str(self) -> &'static str {
        match self {
            AlleleEventClass::Reference => "reference",
            AlleleEventClass::Snp => "snp",
            AlleleEventClass::Mnv => "mnv",
            AlleleEventClass::Insertion => "insertion",
            AlleleEventClass::Deletion => "deletion",
            AlleleEventClass::Delins => "delins",
            AlleleEventClass::ComplexIndel => "complex_indel",
            AlleleEventClass::Symbolic => "symbolic",
        }
    }

    pub fn has_indel_component(self) -> bool {
        matches!(
            self,
            AlleleEventClass::Insertion
                | AlleleEventClass::Deletion
                | AlleleEventClass::Delins
                | AlleleEventClass::ComplexIndel
                | AlleleEventClass::Symbolic
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleComponentKind {
    Snp,
    Insertion,
    Deletion,
    Delins,
    Symbolic,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleComponent {
    pub kind: AlleleComponentKind,
    /// 1-based reference coordinate. Insertions use the anchor base coordinate
    /// after which the inserted sequence occurs.
    pub position: usize,
    pub ref_allele: String,
    pub alt_allele: String,
}

impl AlleleComponent {
    pub fn label(&self) -> String {
        match self.kind {
            AlleleComponentKind::Snp => {
                format!(
                    "SNV:{}:{}>{}",
                    self.position, self.ref_allele, self.alt_allele
                )
            }
            AlleleComponentKind::Insertion => {
                format!("INS:{}:+{}", self.position, self.alt_allele)
            }
            AlleleComponentKind::Deletion => {
                let end = self.position + self.ref_allele.len().saturating_sub(1);
                if self.position == end {
                    format!("DEL:{}:{}", self.position, self.ref_allele)
                } else {
                    format!("DEL:{}-{}:{}", self.position, end, self.ref_allele)
                }
            }
            AlleleComponentKind::Delins => {
                format!(
                    "DELINS:{}:{}>{}",
                    self.position, self.ref_allele, self.alt_allele
                )
            }
            AlleleComponentKind::Symbolic => {
                format!(
                    "SYMBOLIC:{}:{}>{}",
                    self.position, self.ref_allele, self.alt_allele
                )
            }
        }
    }
}

pub fn parse_component_label(label: &str) -> Option<AlleleComponent> {
    let (kind, rest) = label.split_once(':')?;
    let (position_raw, allele_raw) = rest.split_once(':')?;
    let position = position_raw
        .split_once('-')
        .map(|(start, _)| start)
        .unwrap_or(position_raw)
        .parse::<usize>()
        .ok()?;

    match kind {
        "SNV" => {
            let (ref_allele, alt_allele) = allele_raw.split_once('>')?;
            Some(AlleleComponent {
                kind: AlleleComponentKind::Snp,
                position,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
            })
        }
        "INS" => Some(AlleleComponent {
            kind: AlleleComponentKind::Insertion,
            position,
            ref_allele: String::new(),
            alt_allele: allele_raw.strip_prefix('+')?.to_string(),
        }),
        "DEL" => Some(AlleleComponent {
            kind: AlleleComponentKind::Deletion,
            position,
            ref_allele: allele_raw.to_string(),
            alt_allele: String::new(),
        }),
        "DELINS" => {
            let (ref_allele, alt_allele) = allele_raw.split_once('>')?;
            Some(AlleleComponent {
                kind: AlleleComponentKind::Delins,
                position,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
            })
        }
        "SYMBOLIC" => {
            let (ref_allele, alt_allele) = allele_raw.split_once('>')?;
            Some(AlleleComponent {
                kind: AlleleComponentKind::Symbolic,
                position,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
            })
        }
        _ => None,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleEvent {
    pub class: AlleleEventClass,
    pub components: Vec<AlleleComponent>,
    pub affected_start: usize,
    pub affected_end: usize,
}

impl AlleleEvent {
    pub fn component_labels(&self) -> Vec<String> {
        self.components.iter().map(AlleleComponent::label).collect()
    }
}

fn slice_chars(chars: &[char], start: usize, end: usize) -> String {
    chars[start..end].iter().collect()
}

fn classify_equal_length(position: usize, ref_chars: &[char], alt_chars: &[char]) -> AlleleEvent {
    let mut components = Vec::new();
    for (idx, (&ref_base, &alt_base)) in ref_chars.iter().zip(alt_chars.iter()).enumerate() {
        if !ref_base.eq_ignore_ascii_case(&alt_base) {
            components.push(AlleleComponent {
                kind: AlleleComponentKind::Snp,
                position: position + idx,
                ref_allele: ref_base.to_string(),
                alt_allele: alt_base.to_string(),
            });
        }
    }

    let class = match components.len() {
        0 => AlleleEventClass::Reference,
        1 => AlleleEventClass::Snp,
        _ => AlleleEventClass::Mnv,
    };
    AlleleEvent {
        class,
        components,
        affected_start: position,
        affected_end: position + ref_chars.len().saturating_sub(1),
    }
}

/// Decompose a VCF-compatible allele into a coarse event class plus simple
/// components. This is intentionally conservative: when a length-changing
/// allele also carries substitutions, it is classified as `complex_indel`
/// rather than flattening the downstream coding effect into false MNVs.
pub fn decompose_allele(position: usize, ref_allele: &str, alt_allele: &str) -> AlleleEvent {
    let ref_chars: Vec<char> = ref_allele.chars().collect();
    let alt_chars: Vec<char> = alt_allele.chars().collect();
    let affected_end = position + ref_chars.len().saturating_sub(1);

    if alt_allele.starts_with('<') && alt_allele.ends_with('>') {
        return AlleleEvent {
            class: AlleleEventClass::Symbolic,
            components: vec![AlleleComponent {
                kind: AlleleComponentKind::Symbolic,
                position,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
            }],
            affected_start: position,
            affected_end,
        };
    }

    if ref_chars.len() == alt_chars.len() {
        return classify_equal_length(position, &ref_chars, &alt_chars);
    }

    let min_len = ref_chars.len().min(alt_chars.len());
    let mut prefix = 0usize;
    while prefix < min_len && ref_chars[prefix].eq_ignore_ascii_case(&alt_chars[prefix]) {
        prefix += 1;
    }

    let mut ref_end = ref_chars.len();
    let mut alt_end = alt_chars.len();
    while ref_end > prefix
        && alt_end > prefix
        && ref_chars[ref_end - 1].eq_ignore_ascii_case(&alt_chars[alt_end - 1])
    {
        ref_end -= 1;
        alt_end -= 1;
    }

    let event_pos = position + prefix;
    let ref_core = slice_chars(&ref_chars, prefix, ref_end);
    let alt_core = slice_chars(&alt_chars, prefix, alt_end);

    if ref_core.is_empty() && !alt_core.is_empty() {
        return AlleleEvent {
            class: AlleleEventClass::Insertion,
            components: vec![AlleleComponent {
                kind: AlleleComponentKind::Insertion,
                position: event_pos.saturating_sub(1).max(position),
                ref_allele: String::new(),
                alt_allele: alt_core,
            }],
            affected_start: position,
            affected_end,
        };
    }

    if alt_core.is_empty() && !ref_core.is_empty() {
        return AlleleEvent {
            class: AlleleEventClass::Deletion,
            components: vec![AlleleComponent {
                kind: AlleleComponentKind::Deletion,
                position: event_pos,
                ref_allele: ref_core.clone(),
                alt_allele: String::new(),
            }],
            affected_start: event_pos,
            affected_end: event_pos + ref_core.len().saturating_sub(1),
        };
    }

    let ref_core_chars: Vec<char> = ref_core.chars().collect();
    let alt_core_chars: Vec<char> = alt_core.chars().collect();
    let core_min = ref_core_chars.len().min(alt_core_chars.len());
    let mut components = Vec::new();

    for idx in 0..core_min {
        let ref_base = ref_core_chars[idx];
        let alt_base = alt_core_chars[idx];
        if !ref_base.eq_ignore_ascii_case(&alt_base) {
            components.push(AlleleComponent {
                kind: AlleleComponentKind::Snp,
                position: event_pos + idx,
                ref_allele: ref_base.to_string(),
                alt_allele: alt_base.to_string(),
            });
        }
    }

    if ref_core_chars.len() > core_min {
        components.push(AlleleComponent {
            kind: AlleleComponentKind::Deletion,
            position: event_pos + core_min,
            ref_allele: slice_chars(&ref_core_chars, core_min, ref_core_chars.len()),
            alt_allele: String::new(),
        });
    } else if alt_core_chars.len() > core_min {
        components.push(AlleleComponent {
            kind: AlleleComponentKind::Insertion,
            position: event_pos + core_min.saturating_sub(1),
            ref_allele: String::new(),
            alt_allele: slice_chars(&alt_core_chars, core_min, alt_core_chars.len()),
        });
    }

    let class = if components
        .iter()
        .any(|c| matches!(c.kind, AlleleComponentKind::Snp))
        && components.iter().any(|c| {
            matches!(
                c.kind,
                AlleleComponentKind::Insertion | AlleleComponentKind::Deletion
            )
        }) {
        AlleleEventClass::ComplexIndel
    } else {
        AlleleEventClass::Delins
    };

    if components.is_empty() {
        components.push(AlleleComponent {
            kind: AlleleComponentKind::Delins,
            position: event_pos,
            ref_allele: ref_core.clone(),
            alt_allele: alt_core.clone(),
        });
    }

    AlleleEvent {
        class,
        components,
        affected_start: event_pos,
        affected_end: event_pos + ref_core.len().saturating_sub(1),
    }
}

pub fn substitution_components(
    position: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> Vec<AlleleComponent> {
    let event = decompose_allele(position, ref_allele, alt_allele);
    if event.class.has_indel_component() {
        Vec::new()
    } else {
        event.components
    }
}

#[cfg(test)]
mod tests {
    use super::{decompose_allele, parse_component_label, AlleleComponentKind, AlleleEventClass};

    #[test]
    fn decomposes_mnv() {
        let event = decompose_allele(10, "AC", "GT");
        assert_eq!(event.class, AlleleEventClass::Mnv);
        assert_eq!(event.components.len(), 2);
        assert_eq!(event.components[0].position, 10);
        assert_eq!(event.components[1].position, 11);
    }

    #[test]
    fn decomposes_anchored_insertion() {
        let event = decompose_allele(10, "A", "ATG");
        assert_eq!(event.class, AlleleEventClass::Insertion);
        assert_eq!(event.components[0].kind, AlleleComponentKind::Insertion);
        assert_eq!(event.components[0].position, 10);
        assert_eq!(event.components[0].alt_allele, "TG");
    }

    #[test]
    fn decomposes_anchored_deletion() {
        let event = decompose_allele(10, "ATG", "A");
        assert_eq!(event.class, AlleleEventClass::Deletion);
        assert_eq!(event.components[0].kind, AlleleComponentKind::Deletion);
        assert_eq!(event.components[0].position, 11);
        assert_eq!(event.components[0].ref_allele, "TG");
    }

    #[test]
    fn decomposes_complex_indel_with_substitution() {
        let event = decompose_allele(10, "TC", "G");
        assert_eq!(event.class, AlleleEventClass::ComplexIndel);
        assert_eq!(event.components.len(), 2);
        assert_eq!(event.component_labels(), vec!["SNV:10:T>G", "DEL:11:C"]);
    }

    #[test]
    fn parses_component_labels_roundtrip() {
        for label in [
            "SNV:10:A>G",
            "INS:10:+T",
            "DEL:11-12:TG",
            "DELINS:20:AC>G",
            "SYMBOLIC:30:N><DEL>",
        ] {
            let component = parse_component_label(label).expect("component label");
            assert_eq!(component.label(), label);
        }
    }
}
