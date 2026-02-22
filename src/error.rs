use serde::Serialize;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs::File;
use std::io::Write;
use thiserror::Error;

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
        }
    }
}

impl Display for ErrorCode {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Error)]
#[error("[{code}] {message}")]
pub struct AppError {
    pub code: ErrorCode,
    pub message: String,
}

impl AppError {
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
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
}

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
        AppError::new(ErrorCode::Io, format!("I/O error: {}", value))
    }
}

impl From<csv::Error> for AppError {
    fn from(value: csv::Error) -> Self {
        AppError::new(ErrorCode::Csv, format!("CSV error: {}", value))
    }
}

impl From<rust_htslib::errors::Error> for AppError {
    fn from(value: rust_htslib::errors::Error) -> Self {
        AppError::new(ErrorCode::Htslib, format!("HTSlib error: {}", value))
    }
}

impl From<std::str::Utf8Error> for AppError {
    fn from(value: std::str::Utf8Error) -> Self {
        AppError::new(ErrorCode::Utf8, format!("UTF-8 error: {}", value))
    }
}

impl From<std::string::FromUtf8Error> for AppError {
    fn from(value: std::string::FromUtf8Error) -> Self {
        AppError::new(
            ErrorCode::Utf8,
            format!("UTF-8 conversion error: {}", value),
        )
    }
}

impl From<std::num::ParseIntError> for AppError {
    fn from(value: std::num::ParseIntError) -> Self {
        AppError::new(
            ErrorCode::ParseInt,
            format!("Integer parsing error: {}", value),
        )
    }
}

pub type AppResult<T> = Result<T, AppError>;

#[derive(Debug, Serialize)]
struct ErrorJson<'a> {
    schema_version: &'static str,
    code: &'a str,
    exit_code: i32,
    message: &'a str,
}

pub fn error_to_json(error: &AppError) -> String {
    let payload = ErrorJson {
        schema_version: "1.0.0",
        code: error.code.as_str(),
        exit_code: error.code.exit_code(),
        message: &error.message,
    };
    serde_json::to_string_pretty(&payload)
        .unwrap_or_else(|_| "{\"schema_version\":\"1.0.0\",\"code\":\"E000\",\"exit_code\":1,\"message\":\"Failed to serialize error\"}".to_string())
        + "\n"
}

pub fn write_error_json(path: &str, error: &AppError) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    let payload = ErrorJson {
        schema_version: "1.0.0",
        code: error.code.as_str(),
        exit_code: error.code.exit_code(),
        message: &error.message,
    };
    serde_json::to_writer_pretty(&mut file, &payload)?;
    file.write_all(b"\n")
}

#[cfg(test)]
mod tests {
    use super::{AppError, ErrorCode};

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
}
