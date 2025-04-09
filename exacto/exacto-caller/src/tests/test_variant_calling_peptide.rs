// extern crate polars;
//
// use polars::prelude::*;
// use std::fs;
// use std::path::Path;
// use crate::algorithms::variant_calling_peptide::identify_peptide_variants;
// use crate::structs::mutant_peptides_set::MutantPeptidesSet;
//
//
// #[test]
// fn test_identify_peptide_variants_1() {
//     let fasta_file_path = Path::new("src/tests/data/fasta/sample_peptide_sequences.fa.gz");
//     let fasta_full_path = fs::canonicalize(fasta_file_path).unwrap();
//     let fasta_file: &str = fasta_full_path.to_str().unwrap();
//     let reference_fasta_file_path = Path::new("src/tests/data/fasta/reference_peptide_sequences.fa.gz");
//     let reference_fasta_full_path = fs::canonicalize(reference_fasta_file_path).unwrap();
//     let reference_fasta_file: &str = reference_fasta_full_path.to_str().unwrap();
//     let translations_tsv_file_path = Path::new("src/tests/data/tsv/translations/sample_translations.tsv.gz");
//     let translations_tsv_full_path = fs::canonicalize(translations_tsv_file_path).unwrap();
//     let translations_tsv_file: &str = translations_tsv_full_path.to_str().unwrap();
//     let rna_vars_tsv_file_path = Path::new("src/tests/data/tsv/variants/sample_rna_variants.tsv.gz");
//     let rna_vars_tsv_full_path = fs::canonicalize(rna_vars_tsv_file_path).unwrap();
//     let rna_vars_tsv_file: &str = rna_vars_tsv_full_path.to_str().unwrap();
//     let dna_vars_tsv_file_path = Path::new("src/tests/data/tsv/variants/sample_dna_variants.tsv.gz");
//     let dna_vars_tsv_full_path = fs::canonicalize(dna_vars_tsv_file_path).unwrap();
//     let dna_vars_tsv_file: &str = dna_vars_tsv_full_path.to_str().unwrap();
//
//     let mutant_peptides_set: MutantPeptidesSet = identify_peptide_variants(
//         fasta_file,
//         "",
//         "",
//         reference_fasta_file,
//         translations_tsv_file,
//         rna_vars_tsv_file,
//         dna_vars_tsv_file,
//         "",
//         3,
//         8,
//         1,
//         1000
//     );
//
//     assert!(mutant_peptides_set.get_size() == 1);
//     for mutant_peptide in mutant_peptides_set.mutant_peptides.iter() {
//         if mutant_peptide.peptide_sequence == "MCDEFGHIKLMNPQRSTVWY".into() {
//             assert!(mutant_peptide.rna_variants_read_names[0].contains("read1"));
//             assert!(mutant_peptide.rna_variants_read_names[0].contains("read2"));
//             assert!(mutant_peptide.rna_variants_read_names[0].contains("read3"));
//             assert!(mutant_peptide.dna_variants_read_names[0].contains("read100"));
//             assert!(mutant_peptide.dna_variants_read_names[0].contains("read101"));
//             assert!(mutant_peptide.dna_variants_read_names[0].contains("read102"));
//             assert!(mutant_peptide.rna_variants[0] == "chr1:100:+:D:chr1:101:+:U:INS:ACGATGCTAGCTAGTCGATCGTAGC:25".into());
//             assert!(mutant_peptide.dna_variants[0] == "chr1:100:+:D:chr1:101:+:U:INS:ACGATGCTAGCTAGTCGATCGTAGC:25".into());
//         }
//     }
// }
