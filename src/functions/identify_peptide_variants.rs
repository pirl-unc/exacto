// /// This function identifies peptides variants.
// #[pyfunction]
// fn identify_peptide_variants(
//     py: Python,
//     fasta_file: String,
//     rna_bam_file: String,
//     rna_bam_bai_file: String,
//     reference_fasta_file: String,
//     translations_tsv_file: String,
//     rna_variants_tsv_file: String,
//     dna_variants_tsv_file: String,
//     exclude_bed_file: String,
//     min_reads: usize,
//     k: usize,
//     num_threads: usize,
//     dna_variant_padding: usize,
//     output_tsv_file: String,
//     output_fasta_file: String,
//     gzip: bool,
//     output_type: String
// ) -> PyResult<PyDataFrame> {
//     let mutant_peptides_set: caller::MutantPeptidesSet = caller::identify_peptide_variants(
//         fasta_file.as_str(),
//         rna_bam_file.as_str(),
//         rna_bam_bai_file.as_str(),
//         reference_fasta_file.as_str(),
//         translations_tsv_file.as_str(),
//         rna_variants_tsv_file.as_str(),
//         dna_variants_tsv_file.as_str(),
//         exclude_bed_file.as_str(),
//         min_reads,
//         k,
//         num_threads,
//         dna_variant_padding
//     );
//     capture_memory_usage("Successfully ran Exacto peptide variant calling (1)");
//     match output_type.as_str() {
//         "dataframe" => {
//             Ok(PyDataFrame(mutant_peptides_set.to_dataframe(num_threads)))
//         }
//         "file" => {
//             // FASTA file
//             if gzip {
//                 let output_fasta_file_ = File::create(output_fasta_file).unwrap();
//                 let bgzf_writer = bgzf::Writer::new(output_fasta_file_);
//                 let mut fasta_writer = fasta::Writer::new(bgzf_writer);
//                 for mutant_peptide in mutant_peptides_set.mutant_peptides.iter() {
//                     let definition = Definition::new(mutant_peptide.id.to_string(), None);
//                     let sequence = Sequence::from(mutant_peptide.peptide_sequence.as_bytes().to_vec());
//                     let record = fasta::Record::new(definition, sequence);
//                     fasta_writer.write_record(&record).unwrap();
//                 }
//             } else {
//                 let output_fasta_file_ = File::create(output_fasta_file).unwrap();
//                 let mut fasta_writer = fasta::Writer::new(output_fasta_file_);
//                 for mutant_peptide in mutant_peptides_set.mutant_peptides.iter() {
//                     let definition = Definition::new(mutant_peptide.id.to_string(), None);
//                     let sequence = Sequence::from(mutant_peptide.peptide_sequence.as_bytes().to_vec());
//                     let record = fasta::Record::new(definition, sequence);
//                     fasta_writer.write_record(&record).unwrap();
//                 }
//             }
//
//             // TSV file
//             mutant_peptides_set.to_tsv(
//                 output_tsv_file.as_str(),
//                 100_000,
//                 num_threads,
//                 gzip
//             );
//             capture_memory_usage("Successfully wrote to TSV file");
//             Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
//         }
//         other => {
//             let error_message = format!("Unsupported value for output_type: {}", other);
//             Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
//         }
//     }
// }
