//! Utilities for inspecting references stored in MG-index files.

use std::io::Write;

use crate::error::MtsvResult;
use crate::index::MGIndex;
use crate::io::from_file;

/// Write a tab-separated table of the references stored in a list of indices.
///
/// Index numbers are one-based and correspond to the order of `index_paths`.
pub fn write_reference_list<W: Write>(index_paths: &[&str], writer: &mut W) -> MtsvResult<()> {
    writeln!(writer, "index\ttaxid\tgenome_id\tlength")?;

    for (index_number, path) in index_paths.iter().enumerate() {
        info!(
            "Reading reference metadata from index {}: {}",
            index_number + 1,
            path
        );
        let index = from_file::<MGIndex>(path)?;
        for (taxid, genome_id, length) in index.reference_metadata() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                index_number + 1,
                taxid,
                genome_id,
                length
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_reference_list;
    use crate::index::{Database, Gi, MGIndex, TaxId};
    use crate::io::write_to_file;
    use tempfile::NamedTempFile;

    #[test]
    fn writes_metadata_for_multiple_indices() {
        let mut first_db = Database::new();
        first_db.insert(TaxId(9), vec![(Gi(101), b"ACGT".to_vec())]);
        let mut second_db = Database::new();
        second_db.insert(
            TaxId(12),
            vec![(Gi(202), b"AAAAAA".to_vec()), (Gi(203), b"TT".to_vec())],
        );

        let first_file = NamedTempFile::new().unwrap();
        let second_file = NamedTempFile::new().unwrap();
        write_to_file(
            &MGIndex::new(first_db, 2, 2),
            first_file.path().to_str().unwrap(),
        )
        .unwrap();
        write_to_file(
            &MGIndex::new(second_db, 2, 2),
            second_file.path().to_str().unwrap(),
        )
        .unwrap();

        let paths = [
            first_file.path().to_str().unwrap(),
            second_file.path().to_str().unwrap(),
        ];
        let mut output = Vec::new();
        write_reference_list(&paths, &mut output).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "index\ttaxid\tgenome_id\tlength\n1\t9\t101\t4\n2\t12\t202\t6\n2\t12\t203\t2\n"
        );
    }
}
