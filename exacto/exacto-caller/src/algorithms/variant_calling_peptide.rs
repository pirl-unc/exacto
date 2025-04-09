// // Licensed under the Apache License, Version 2.0 (the "License");
// // you may not use this file except in compliance with the License.
// // You may obtain a copy of the License at
// //
// //      http://www.apache.org/licenses/LICENSE-2.0
// //
// // Unless required by applicable law or agreed to in writing, software
// // distributed under the License is distributed on an "AS IS" BASIS,
// // WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// // See the License for the specific language governing permissions and
// // limitations under the License.
//
//
// use bio::io::bed;
// use exacto_util::prelude::*;
// use noodles_bam::{self as bam};
// use polars::prelude::*;
// use rayon::prelude::*;
// use std::collections::{BTreeMap, HashMap, HashSet};
// use std::str;
// use bimap::BiMap;
//
// use crate::common::bam::{create_read_names_map, fetch_bam_records};
// use crate::log_info;
// use crate::structs::mutant_peptide::MutantPeptide;
// use crate::structs::mutant_peptides_set::MutantPeptidesSet;
//
//
// fn unpack_dna_variants(
//     df: &mut DataFrame,
//     bin_size: u32
// ) -> (HashMap<(Box<str>,u32),BTreeMap<u32,Box<str>>>,
//       HashMap<(Box<str>,u32),BTreeMap<u32,Box<str>>>,
//       HashMap<Box<str>,(Box<str>,Box<str>)>) {
//     let df_: DataFrame = df
//         .with_column(
//             df
//                 .column("variant_call_id")
//                 .unwrap()
//                 .cast(&DataType::String).unwrap()
//         )
//         .unwrap()
//         .clone();
//     let mut pos_1_tree_map: HashMap<(Box<str>,u32),BTreeMap<u32,Box<str>>> = HashMap::new();
//     let mut pos_2_tree_map: HashMap<(Box<str>,u32),BTreeMap<u32,Box<str>>> = HashMap::new();
//     let mut variants_map: HashMap<Box<str>,(Box<str>,Box<str>)> = HashMap::new();
//     let row_count = df_.height();
//
//     let variant_call_id_col = df_.column("variant_call_id").unwrap().str().unwrap();
//     let chromosome_1_col = df_.column("chromosome_1").unwrap().str().unwrap();
//     let position_1_col = df_.column("position_1").unwrap().i64().unwrap();
//     let strand_1_col = df_.column("strand_1").unwrap().str().unwrap();
//     let orientation_1_col = df_.column("orientation_1").unwrap().str().unwrap();
//     let chromosome_2_col = df_.column("chromosome_2").unwrap().str().unwrap();
//     let position_2_col = df_.column("position_2").unwrap().i64().unwrap();
//     let strand_2_col = df_.column("strand_2").unwrap().str().unwrap();
//     let orientation_2_col = df_.column("orientation_2").unwrap().str().unwrap();
//     let variant_type_col = df_.column("variant_type").unwrap().str().unwrap();
//     let variant_sequence_col = df_.column("variant_sequence").unwrap().str().unwrap();
//     let consensus_read_names_col = df_.column("consensus_read_names").unwrap().str().unwrap();
//
//     for row_idx in 0..row_count {
//         let mut variant_sequence: Box<str> = "".into();
//         let mut variant_sequence_length: usize = 0;
//         if variant_sequence_col.get(row_idx).is_some() {
//             variant_sequence = variant_sequence_col.get(row_idx).unwrap().to_string().into_boxed_str();
//             variant_sequence_length = variant_sequence.len();
//         }
//
//         let variant = format!(
//             "{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
//             chromosome_1_col.get(row_idx).unwrap(),
//             position_1_col.get(row_idx).unwrap(),
//             strand_1_col.get(row_idx).unwrap(),
//             orientation_1_col.get(row_idx).unwrap(),
//             chromosome_2_col.get(row_idx).unwrap(),
//             position_2_col.get(row_idx).unwrap(),
//             strand_2_col.get(row_idx).unwrap(),
//             orientation_2_col.get(row_idx).unwrap(),
//             variant_type_col.get(row_idx).unwrap(),
//             variant_sequence,
//             variant_sequence_length
//         );
//         let variant_call_id: Box<str> = variant_call_id_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let read_names: Box<str> = consensus_read_names_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let chromosome_1: Box<str> = chromosome_1_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let chromosome_2: Box<str> = chromosome_2_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let position_1: u32 = position_1_col.get(row_idx).unwrap() as u32;
//         let position_2: u32 = position_2_col.get(row_idx).unwrap() as u32;
//         let zipcode_1: u32 = position_1 / bin_size;
//         let zipcode_2: u32 = position_2 / bin_size;
//
//         pos_1_tree_map
//             .entry((chromosome_1.clone(),zipcode_1))
//             .or_insert_with(BTreeMap::new)
//             .insert(position_1,variant_call_id.clone());
//         pos_2_tree_map
//             .entry((chromosome_2.clone(),zipcode_2))
//             .or_insert_with(BTreeMap::new)
//             .insert(position_2,variant_call_id.clone());
//         variants_map.insert(variant_call_id, (variant.into_boxed_str(),read_names));
//     }
//
//     (pos_1_tree_map,pos_2_tree_map,variants_map)
// }
//
// fn unpack_rna_variants(
//     df_rna_variants: &mut DataFrame
// ) -> (HashMap<Box<str>,Vec<(Box<str>,Box<str>)>>,
//       HashMap<Box<str>,(Box<str>,u32,Box<str>,u32)>) {
//     let df_rna_variants_: DataFrame = df_rna_variants
//         .with_column(
//             df_rna_variants
//                 .column("variant_call_id")
//                 .unwrap()
//                 .cast(&DataType::String).unwrap()
//         )
//         .unwrap()
//         .clone();
//
//     let mut rna_variant_call_ids_map: HashMap<Box<str>,(Box<str>,u32,Box<str>,u32)> = HashMap::new();
//     let mut rna_read_names_variants_map: HashMap<Box<str>,Vec<(Box<str>,Box<str>)>> = HashMap::new();
//     let row_count = df_rna_variants.height();
//
//     let variant_call_id_col = df_rna_variants_.column("variant_call_id").unwrap().str().unwrap();
//     let chromosome_1_col = df_rna_variants_.column("chromosome_1").unwrap().str().unwrap();
//     let position_1_col = df_rna_variants_.column("position_1").unwrap().i64().unwrap();
//     let strand_1_col = df_rna_variants_.column("strand_1").unwrap().str().unwrap();
//     let orientation_1_col = df_rna_variants_.column("orientation_1").unwrap().str().unwrap();
//     let chromosome_2_col = df_rna_variants_.column("chromosome_2").unwrap().str().unwrap();
//     let position_2_col = df_rna_variants_.column("position_2").unwrap().i64().unwrap();
//     let strand_2_col = df_rna_variants_.column("strand_2").unwrap().str().unwrap();
//     let orientation_2_col = df_rna_variants_.column("orientation_2").unwrap().str().unwrap();
//     let variant_type_col = df_rna_variants_.column("variant_type").unwrap().str().unwrap();
//     let variant_sequence_col = df_rna_variants_.column("variant_sequence").unwrap().str().unwrap();
//     let consensus_read_names_col = df_rna_variants_.column("consensus_read_names").unwrap().str().unwrap();
//
//     for row_idx in 0..row_count {
//         let mut variant_sequence: Box<str> = "".into();
//         let mut variant_sequence_length: usize = 0;
//         if variant_sequence_col.get(row_idx).is_some() {
//             variant_sequence = variant_sequence_col.get(row_idx).unwrap().to_string().into_boxed_str();
//             variant_sequence_length = variant_sequence.len();
//         }
//
//         let rna_variant = format!(
//             "{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
//             chromosome_1_col.get(row_idx).unwrap(),
//             position_1_col.get(row_idx).unwrap(),
//             strand_1_col.get(row_idx).unwrap(),
//             orientation_1_col.get(row_idx).unwrap(),
//             chromosome_2_col.get(row_idx).unwrap(),
//             position_2_col.get(row_idx).unwrap(),
//             strand_2_col.get(row_idx).unwrap(),
//             orientation_2_col.get(row_idx).unwrap(),
//             variant_type_col.get(row_idx).unwrap(),
//             variant_sequence,
//             variant_sequence_length
//         );
//
//         let variant_call_id: Box<str> = variant_call_id_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let chromosome_1: Box<str> = chromosome_1_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let chromosome_2: Box<str> = chromosome_2_col.get(row_idx).unwrap().to_string().into_boxed_str();
//         let position_1: u32 = position_1_col.get(row_idx).unwrap() as u32;
//         let position_2: u32 = position_2_col.get(row_idx).unwrap() as u32;
//
//         rna_variant_call_ids_map.insert(variant_call_id.clone(), (chromosome_1,position_1,chromosome_2,position_2));
//
//         if let Some(consensus_read_names) = consensus_read_names_col.get(row_idx) {
//             for read_name in consensus_read_names.split(',') {
//                 rna_read_names_variants_map
//                     .entry(read_name.to_string().into_boxed_str())
//                     .or_default()
//                     .push((variant_call_id.clone(), rna_variant.to_string().into_boxed_str()));
//             }
//         }
//     }
//
//     (rna_read_names_variants_map,rna_variant_call_ids_map)
// }
//
// /// Identify mutant k-mer sequences.
// ///
// /// Returns:
// /// * HashMap<peptide_id,(peptide_sequence,kmer_mask,rna_read_names)>
// fn identify_mutant_kmer_sequences(
//     fasta_file: &str,
//     reference_fasta_file: &str,
//     translations_tsv_file: &str,
//     rna_variant_read_names: HashSet<Box<str>>,
//     min_reads: usize,
//     k: usize,
//     num_threads: usize
// ) -> HashMap<Box<str>,(Box<str>,Vec<bool>,HashSet<Box<str>>)> {
//     let thread_pool = rayon::ThreadPoolBuilder::new()
//         .num_threads(num_threads)
//         .build()
//         .unwrap();
//
//     // Step 1: Read translations and filter peptides with enough reads
//     let mut df_translations: DataFrame = read_tsv_file(translations_tsv_file);
//     df_translations = df_translations
//         .with_column(
//             df_translations
//                 .column("peptide_id")
//                 .unwrap()
//                 .cast(&DataType::String).unwrap()
//             )
//         .unwrap()
//         .clone();
//     let df_read_names = DataFrame::new(vec![
//         Column::from(Series::new("read_name".into(), rna_variant_read_names.iter().map(|s| s.as_ref()).collect::<Vec<_>>()))
//     ]).unwrap();
//     let supported_peptide_ids: HashSet<Box<str>> = df_translations
//         .clone()
//         .lazy()
//         .join(
//             df_read_names.lazy(),
//             [col("read_name")],
//             [col("read_name")],
//             JoinArgs::new(JoinType::Inner)
//         )
//         .group_by([col("peptide_id")])
//         .agg([col("read_name").n_unique().alias("unique_read_count")])
//         .filter(col("unique_read_count").gt_eq(lit(min_reads as i64)))
//         .collect()
//         .unwrap()
//         .column("peptide_id")
//         .unwrap()
//         .str()
//         .unwrap()
//         .into_iter()
//         .filter_map(|v| v.map(|id| id.to_string().into_boxed_str()))
//         .collect();
//
//     // Step 2: Load and filter peptide sequences
//     let mut peptide_sequences: HashMap<Box<str>,Box<str>> = read_fasta_file(fasta_file).into_iter().collect();
//     peptide_sequences.retain(|key, _| supported_peptide_ids.contains(key));
//     let reference_peptide_sequences: Vec<(Box<str>,Box<str>)> = read_fasta_file(reference_fasta_file);
//
//     // Step 3. Identify k-mers
//     let query_kmers: HashSet<Box<str>> = thread_pool.install(|| {
//         peptide_sequences
//             .par_iter()
//             .flat_map(|(id,sequence)| {
//                 if sequence.len() >= k {
//                     find_kmers(&**sequence, k)
//                         .into_keys()
//                         .collect::<Vec<_>>()
//                 } else {
//                     Vec::new()
//                 }
//             })
//             .collect()
//     });
//     let reference_kmers: HashSet<Box<str>> = thread_pool.install(|| {
//         reference_peptide_sequences
//             .par_iter()
//             .flat_map(|(id,sequence)| {
//                 if sequence.len() >= k {
//                     find_kmers(&**sequence, k)
//                         .into_keys()
//                         .collect::<Vec<_>>()
//                 } else {
//                     Vec::new()
//                 }
//             })
//             .collect()
//     });
//     let unique_kmers: HashSet<Box<str>> = query_kmers.difference(&reference_kmers).cloned().collect();
//
//     // Step 4: Identify mutant peptides
//     let mutant_peptides: Vec<(Box<str>, Box<str>, Vec<bool>)> = thread_pool.install(|| {
//         peptide_sequences
//             .par_iter()
//             .filter_map(|(peptide_id, sequence)| {
//                 if sequence.len() >= k {
//                     let kmer_mask: Vec<bool> = (0..=sequence.len() - k)
//                         .map(|i| {
//                             let kmer: Box<str> = sequence[i..i + k].to_string().into_boxed_str();
//                             unique_kmers.contains(&kmer)
//                         })
//                         .collect();
//                     if kmer_mask.iter().any(|&is_unique| is_unique) {
//                         Some((peptide_id.clone(), sequence.clone(), kmer_mask))
//                     } else {
//                         None
//                     }
//                 } else {
//                     None
//                 }
//             })
//             .collect()
//     });
//
//     // Step 5: Create mutant peptide map
//     mutant_peptides
//         .into_iter()
//         .map(|(peptide_id, sequence, kmer_mask)| {
//             let rna_read_names: HashSet<Box<str>> = df_translations
//                 .clone()
//                 .lazy()
//                 .filter(col("peptide_id").eq(lit(&*peptide_id)))
//                 .collect()
//                 .unwrap()
//                 .column("read_name")
//                 .unwrap()
//                 .str()
//                 .unwrap()
//                 .into_iter()
//                 .filter_map(|v| v.map(|s| s.to_string().into_boxed_str()))
//                 .collect();
//             (peptide_id, (sequence, kmer_mask, rna_read_names))
//         })
//         .collect()
// }
//
// pub fn identify_peptide_variants(
//     fasta_file: &str,
//     rna_bam_file: &str,
//     rna_bam_bai_file: &str,
//     reference_fasta_file: &str,
//     translations_tsv_file: &str,
//     rna_variants_tsv_file: &str,
//     dna_variants_tsv_file: &str,
//     exclude_bed_file: &str,
//     min_reads: usize,
//     k: usize,
//     num_threads: usize,
//     dna_variant_padding: usize
// ) -> MutantPeptidesSet {
//     let dna_variant_padding_: u32 = dna_variant_padding as u32;
//
//     // Step 1. Identify read names to exclude based on excluded BED file
//     let mut excluded_rna_read_names: HashSet<Box<str>> = HashSet::new();
//     if exclude_bed_file.is_empty() == false {
//         let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
//             rna_bam_file,
//             rna_bam_bai_file,
//             num_threads
//         );
//         let bed_records: Vec<bed::Record> = read_bed_file(exclude_bed_file);
//         for bed_record in bed_records.iter() {
//             let mut records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
//                 rna_bam_file,
//                 rna_bam_bai_file,
//                 bed_record.chrom(),
//                 bed_record.start() as usize,
//                 bed_record.end() as usize,
//                 &read_names_map,
//                 num_threads
//             );
//             for read_id in records_map.keys() {
//                 let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();
//                 excluded_rna_read_names.insert(read_name);
//             }
//         }
//     }
//
//     // Step 2. Load RNA and DNA variants
//     log_info!("Reading the DNA variants");
//     let mut df_rna_variants: DataFrame = read_tsv_file(rna_variants_tsv_file);
//     log_info!("Reading the RNA variants");
//     let mut df_dna_variants: DataFrame = read_tsv_file(dna_variants_tsv_file);
//
//     // Step 3. Unpack RNA variants
//     log_info!("Unpacking the RNA variants");
//     let (rna_read_names_variants_map,rna_variant_call_ids_map) = unpack_rna_variants(&mut df_rna_variants);
//
//     // Step 4. Build binary trees of the DNA variants
//     log_info!("Unpacking the DNA variants");
//     let bin_size: u32 = 10000;
//     let (dna_pos_1_map, dna_pos_2_map, dna_variants_map) = unpack_dna_variants(&mut df_dna_variants, bin_size);
//
//     // Step 5. Identify mutant k-mer sequences
//     log_info!("Identifying mutant {}-mer sequences", k);
//     let mut mutant_peptides_map: HashMap<Box<str>,(Box<str>,Vec<bool>,HashSet<Box<str>>)> = identify_mutant_kmer_sequences(
//         fasta_file,
//         reference_fasta_file,
//         translations_tsv_file,
//         rna_read_names_variants_map.keys().cloned().collect(),
//         min_reads,
//         k,
//         num_threads
//     );
//     log_info!("{} mutant peptide sequences with at least 1 unique {}-mers", mutant_peptides_map.len(), k);
//
//     // Step 6. Identify mutant peptides
//     log_info!("Identifying mutant peptides");
//     let thread_pool = rayon::ThreadPoolBuilder::new()
//         .num_threads(num_threads)
//         .build()
//         .unwrap();
//     let mutant_peptides: Vec<MutantPeptide> = thread_pool.install(|| {
//         mutant_peptides_map
//             .par_iter()
//             .filter_map(|(peptide_id, (peptide_sequence, kmer_mask, rna_read_names))| {
//                 // Make sure none of the RNA read names is in the excluded RNA read names
//                 for rna_read_name in rna_read_names.iter() {
//                     if excluded_rna_read_names.contains(rna_read_name) {
//                         return None;
//                     }
//                 }
//
//                 // Get RNA variants
//                 let mut rna_variant_call_ids: HashSet<Box<str>> = HashSet::new();
//                 let mut rna_variants_read_names_map: HashMap<(Box<str>,Box<str>), usize> = HashMap::new();
//                 for rna_read_name in rna_read_names.iter() {
//                     if let Some(variants) = rna_read_names_variants_map.get(rna_read_name) {
//                         for (variant_call_id, rna_variant) in variants {
//                             *rna_variants_read_names_map
//                                 .entry((variant_call_id.clone(), rna_variant.clone()))
//                                 .or_insert(0) += 1;
//                         }
//                     }
//                 }
//                 let mut rna_variants = Vec::new();
//                 let mut rna_variants_consensus_read_names = Vec::new();
//                 for ((variant_call_id, rna_variant), count) in rna_variants_read_names_map.iter() {
//                     if *count == rna_read_names.len() {
//                         rna_variants.push(rna_variant.clone());
//                         rna_variants_consensus_read_names.push(
//                             rna_read_names
//                                 .iter()
//                                 .map(|s| s.as_ref())
//                                 .collect::<Vec<_>>()
//                                 .join(",")
//                                 .into_boxed_str(),
//                         );
//                         rna_variant_call_ids.insert(variant_call_id.clone());
//                     }
//                 }
//
//                 // Get DNA variants
//                 let mut dna_variants = Vec::new();
//                 let mut dna_variants_consensus_read_names = Vec::new();
//                 for rna_variant_call_id in rna_variant_call_ids.iter() {
//                     let (chromosome_1,position_1,chromosome_2,position_2) = rna_variant_call_ids_map.get(rna_variant_call_id).unwrap();
//                     let zipcode_1: u32 = *position_1 / bin_size;
//                     let zipcode_2: u32 = *position_2 / bin_size;
//                     let min_pos_1: u32 = if *position_1 >= dna_variant_padding_ { *position_1 - dna_variant_padding_ } else { 0 };
//                     let max_pos_1: u32 = *position_1 + dna_variant_padding_;
//                     let min_pos_2: u32 = if *position_2 >= dna_variant_padding_ { *position_2 - dna_variant_padding_ } else { 0 };
//                     let max_pos_2: u32 = *position_2 + dna_variant_padding_;
//                     let mut dna_variant_call_ids: HashSet<Box<str>> = HashSet::new();
//                     if dna_pos_1_map.contains_key(&(chromosome_1.to_string().into(),zipcode_1)) {
//                         for (_,dna_variant_call_id) in dna_pos_1_map.get(&(chromosome_1.to_string().into(),zipcode_1)).unwrap().range(min_pos_1..=max_pos_1) {
//                             dna_variant_call_ids.insert(dna_variant_call_id.clone());
//                         }
//                     }
//                     if dna_pos_2_map.contains_key(&(chromosome_2.to_string().into(),zipcode_2)) {
//                         for (_,dna_variant_call_id) in dna_pos_2_map.get(&(chromosome_2.to_string().into(),zipcode_2)).unwrap().range(min_pos_2..=max_pos_2) {
//                             dna_variant_call_ids.insert(dna_variant_call_id.clone());
//                         }
//                     }
//                     for dna_variant_call_id in dna_variant_call_ids.iter() {
//                         dna_variants.push(dna_variants_map.get(dna_variant_call_id).unwrap().0.clone());
//                         dna_variants_consensus_read_names.push(dna_variants_map.get(dna_variant_call_id).unwrap().1.clone());
//                     }
//                 }
//                 if !rna_variants.is_empty() && !dna_variants.is_empty() {
//                     Some(MutantPeptide::new(
//                         peptide_id.clone(),
//                         peptide_sequence.clone(),
//                         rna_read_names.iter().cloned().collect(),
//                         k,
//                         kmer_mask.clone(),
//                         rna_variants,
//                         rna_variants_consensus_read_names,
//                         dna_variants,
//                         dna_variants_consensus_read_names,
//                     ))
//                 } else {
//                     None
//                 }
//             })
//             .collect::<Vec<MutantPeptide>>()
//     });
//
//     log_info!("Preparing a mutant peptides set");
//     let mut mutant_peptides_set: MutantPeptidesSet = MutantPeptidesSet::new();
//     for mutant_peptide in mutant_peptides {
//         mutant_peptides_set.add_mutant_peptide(mutant_peptide);
//     }
//     mutant_peptides_set
// }
