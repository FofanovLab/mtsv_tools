//! Miscellaneous.

use crate::error::*;
use crate::index::{Gi, TaxId};
use chrono::Local;
use env_logger::LogBuilder;
use log::{LogLevelFilter, LogRecord};

/// Initialize the program-wide logger to write to stdout with timestamps.
pub fn init_logging(level: LogLevelFilter) {
    let mut builder = LogBuilder::new();

    builder.filter(None, level).format(|record: &LogRecord| {
        format!(
            "[{} {} {}] {}",
            record.level(),
            Local::now().format("%Y-%m-%d %H:%M:%S%.3f"),
            record.location().module_path(),
            record.args()
        )
    });

    let _ = builder.init();
}

/// Parse a reference sequence's read header in the format expected by mtsv: `ACCESSION-TAXID`.
pub fn parse_read_header(h: &str) -> MtsvResult<(Gi, TaxId)> {
    let (gi, tax_id, alternate) = parse_read_header_with_alternate(h)?;
    if alternate.is_some() {
        return Err(MtsvError::InvalidHeader(String::from(h)));
    }
    Ok((gi, tax_id))
}

/// Parse `ACCESSION-TAXID` or `ACCESSION-TAXID-ALTERNATE_TAXID`.
pub fn parse_read_header_with_alternate(h: &str) -> MtsvResult<(Gi, TaxId, Option<TaxId>)> {
    let mut tokens = h.split('-');

    let gi = match tokens.next() {
        Some(t) => match t.parse::<Gi>() {
            Ok(t) => t,
            Err(_) => return Err(MtsvError::InvalidInteger(t.to_owned())),
        },
        None => return Err(MtsvError::InvalidHeader(String::from(h))),
    };

    let tax_id = match tokens.next() {
        Some(t) => match t.parse::<TaxId>() {
            Ok(t) => t,
            Err(_) => return Err(MtsvError::InvalidInteger(t.to_owned())),
        },
        None => return Err(MtsvError::InvalidHeader(String::from(h))),
    };

    let alternate = match tokens.next() {
        Some(t) => Some(
            t.parse::<TaxId>()
                .map_err(|_| MtsvError::InvalidInteger(t.to_owned()))?,
        ),
        None => None,
    };
    if tokens.next().is_some() {
        return Err(MtsvError::InvalidHeader(String::from(h)));
    }
    Ok((gi, tax_id, alternate))
}

#[cfg(test)]
mod test {
    use crate::index::{Gi, TaxId};

    use super::{init_logging, parse_read_header, parse_read_header_with_alternate};
    use log::LogLevelFilter;

    #[test]
    fn lines_for_the_line_throne() {
        init_logging(LogLevelFilter::Debug);
    }

    #[test]
    fn success() {
        let (found_gi, found_tax) = parse_read_header("12345-908").unwrap();

        assert_eq!(found_gi, Gi(12345));
        assert_eq!(found_tax, TaxId(908));
    }

    #[test]
    fn success_leading_zeros() {
        let (found_gi, found_tax) = parse_read_header("0001-0002").unwrap();

        assert_eq!(found_gi, Gi(1));
        assert_eq!(found_tax, TaxId(2));
    }

    #[test]
    fn success_with_alternate_taxid() {
        assert_eq!(
            parse_read_header_with_alternate("12345-908-701").unwrap(),
            (Gi(12345), TaxId(908), Some(TaxId(701)))
        );
    }

    #[test]
    #[should_panic]
    fn fail_empty_nodash() {
        let _ = parse_read_header("").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_empty() {
        let _ = parse_read_header("-").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_decimal_gi() {
        let _ = parse_read_header("1.0-543").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_decimal_taxid() {
        let _ = parse_read_header("654981-1.071").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_extra() {
        let _ = parse_read_header("1-2-3").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_non_numeric_gi() {
        let _ = parse_read_header("abc-123").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_non_numeric_taxid() {
        let _ = parse_read_header("123-abc").unwrap();
    }
}
