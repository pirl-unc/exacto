use exacto_core::prelude::*;
use std::fs;
use std::path::Path;

use super::*;


/// `SiteSpan::is_covered_by` is the coverage rule `genotype_rna_variants` applies to every
/// read that is neither a member of a call nor a sibling-insertion carrier, so it must agree,
/// exon by exon, with that rule spelled out as plain interval arithmetic: a substitution or
/// deletion is covered when one exon contains its whole reference span; an insertion is
/// covered when one exon reaches its left flank, an exon ending up to the insert length
/// short of the flank included. Inputs are built to hit the edge cases: several
/// chromosomes interleaved, reads whose exons overlap each other (supplementary alignments
/// do this), insert lengths longer than the flank coordinate itself, sites on a chromosome
/// with no exons at all, and reads with no exons. Configurations come from a fixed-seed LCG
/// so the test is deterministic while still adversarial.
#[test]
fn is_covered_by_matches_the_exon_containment_rule() {
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut next = move |bound: u64| -> u64 {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) % bound
    };

    let num_reads: usize = 150;
    let mut read_exons: Vec<Vec<(u16, u32, u32)>> = Vec::with_capacity(num_reads);
    for _ in 0..num_reads {
        let mut exons: Vec<(u16, u32, u32)> = Vec::new();
        // ~5% of reads have no exons (a read whose model was lost mid-pipeline).
        if next(20) != 0 {
            let num_exons: u64 = 1 + next(4);
            for _ in 0..num_exons {
                let chromosome_id: u16 = next(3) as u16;
                let start: u32 = next(1_000) as u32;
                let length: u32 = 1 + next(400) as u32;
                exons.push((chromosome_id, start, start + length));
            }
        }
        read_exons.push(exons);
    }

    let mut sites: Vec<SiteSpan> = Vec::new();
    for _ in 0..120 {
        // Chromosome 3 exists only as sites, never as exons.
        let chromosome_id: u16 = next(4) as u16;
        let reference_start: u32 = next(1_400) as u32;
        let site: SiteSpan = if next(2) == 0 {
            SiteSpan {
                chromosome_id,
                reference_start,
                reference_end: reference_start + next(60) as u32,
                insertion_length: None,
                side_2: None
            }
        } else {
            SiteSpan {
                chromosome_id,
                reference_start,
                reference_end: reference_start + 1,
                // Occasionally longer than reference_start, so any exon reaching the flank qualifies.
                insertion_length: Some(next(2_000) as u32),
                side_2: None
            }
        };
        sites.push(site);
    }

    let mut num_covered: usize = 0;
    for site in sites.iter() {
        for exons in read_exons.iter() {
            // The rule restated in u64, so the insertion window is plain arithmetic.
            let expected: bool = exons.iter().any(|&(chromosome_id, start, end)| {
                chromosome_id == site.chromosome_id
                    && start <= site.reference_start
                    && match site.insertion_length {
                    Some(insertion_length) => site.reference_start as u64 <= end as u64 + insertion_length as u64,
                    None => site.reference_end <= end
                }
            });
            assert_eq!(
                site.is_covered_by(exons),
                expected,
                "site chr{}:{}-{} insertion_length={:?} against exons {:?}",
                site.chromosome_id,
                site.reference_start,
                site.reference_end,
                site.insertion_length,
                exons
            );
            num_covered += expected as usize;
        }
    }

    // The inputs must exercise both outcomes, or the comparison proves nothing.
    assert!(num_covered > 0);
    assert!(num_covered < num_reads * sites.len());

    // A read with no exons never covers anything.
    let no_exons: Vec<(u16, u32, u32)> = Vec::new();
    for site in sites.iter() {
        assert!(!site.is_covered_by(&no_exons));
    }
}

/// `genotype_rna_variants` must place every cluster read on its allele ladder, in order:
/// a member of the call is Alternate; a non-member carrying a different insertion at the
/// call's left flank is NotCovered, never Reference; a read with an exon spanning the site
/// is Reference; and everything else, spliced across the site, ended before it, on another
/// chromosome, or without a model at all, is NotCovered. scga-mini-rna-002 supplies the models
/// (TP53, all reads aligned to the reverse strand): one read carrying the 12-bp insertion
/// CCCATCCGCCTG at chr17:7,674,224-7,674,225, inside its exon 7,674,181-7,674,290, and one
/// reference read with no event there. Both reads have the same 11 exons, from 7,668,421 to
/// 7,687,490, and no soft clip. The calls are synthetic, so each rung of the ladder can be
/// aimed at a known exon, intron, or transcript end.
#[test]
fn scga_mini_rna_002_genotype_rna_variants_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-002-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );

    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file,
        reference_genome_fasta_file,
        &gene_annotator,
        // None,
        &options,
        1,
        ""
    );

    // The two chosen reads, by id, plus a third id the cluster names but no model backs. Every
    // one of the 399 reads in the BAM is modelled.
    let variant_read_id: usize = *transcript_model_set.read_names_map.get_by_left("scga-mini-rna-002-tumor_chunk_0000/1/ccs").unwrap();
    let reference_read_id: usize = *transcript_model_set.read_names_map.get_by_left("scga-mini-rna-002-tumor_chunk_0000/204/ccs").unwrap();
    let unmodelled_read_id: usize = 1_000_000;
    assert_eq!(transcript_model_set.transcript_models.len(), 399);
    let transcript_models_map: HashMap<usize, Arc<TranscriptModel>> = transcript_model_set
        .transcript_models
        .into_iter()
        .map(|transcript_model| (transcript_model.get_read_id(), Arc::new(transcript_model)))
        .collect();
    assert!(!transcript_models_map.contains_key(&unmodelled_read_id));
    let read_ids: HashSet<usize> = HashSet::from([variant_read_id, reference_read_id, unmodelled_read_id]);

    // Both reads' exons, in genome order.
    let exons: Vec<(u16, u32, u32)> = vec![
        (0, 7668421, 7669690), (0, 7670609, 7670715), (0, 7673535, 7673608), (0, 7673701, 7673837),
        (0, 7674181, 7674290), (0, 7674859, 7674971), (0, 7675053, 7675236), (0, 7675994, 7676272),
        (0, 7676382, 7676403), (0, 7676521, 7676622), (0, 7687377, 7687490)
    ];
    for read_id in [variant_read_id, reference_read_id] {
        let mut read_exons: Vec<(u16, u32, u32)> = transcript_models_map[&read_id]
            .get_exons()
            .iter()
            .map(|exon| (exon.reference_chromosome_id, exon.reference_start, exon.reference_end))
            .collect();
        read_exons.sort();
        assert_eq!(read_exons, exons);
    }

    // The variant read's insertion record, as the ground truth spells it: in read order, so on
    // the reverse strand the reverse complement of the BAM's CAGGCGGATGGG. Its other records are
    // sequencing errors elsewhere.
    let insertion_record: VariantRecord = transcript_models_map[&variant_read_id]
        .get_variant_records()
        .iter()
        .find(|record| *record.get_variant_type() == VariantType::Insertion && record.get_position_1() == 7674224)
        .unwrap()
        .clone();
    assert_eq!(&*insertion_record.get_graph_operation().as_boxed_str(), "0:7674224:-:D:0:7674225:-:U:CCCATCCGCCTG:12:INS");

    // Calls are built from records of reads outside the cluster (ids 99 and 100), so that
    // membership is decided by the one record of the variant read included on purpose.
    let call = |id: usize, operation: GraphOperation, member: Option<VariantRecord>| -> VariantCall {
        let mut variant_records: HashSet<VariantRecord> = HashSet::from([
            VariantRecord::new(99, 10, 10, operation.clone()),
            VariantRecord::new(100, 10, 10, operation)
        ]);
        if let Some(member) = member {
            variant_records.insert(member);
        }
        VariantCall::from_variant_records(id, variant_records, 0, 4, 6, 2)
    };
    let operation = |chromosome: u16, position_1: u32, position_2: u32, sequence: &str, variant_type: VariantType| -> GraphOperation {
        GraphOperation::new(
            chromosome,
            position_1,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome,
            position_2,
            Strand::Forward,
            GraphOperationType::Upstream,
            sequence.into(),
            variant_type
        )
    };
    let variant_calls: Vec<VariantCall> = vec![
        // 1. The insertion itself, with the variant read as a member.
        call(1, insertion_record.get_graph_operation().clone(), Some(insertion_record.clone())),
        // 2. A sibling spelling at the same left flank; nobody in the cluster is a member.
        call(2, operation(0, 7674224, 7674225, "GGGGG", VariantType::Insertion), None),
        // 3. An SNV inside exon 7,674,181-7,674,290 of both reads.
        call(3, operation(0, 7674259, 7674261, "A", VariantType::SingleNucleotideVariant), None),
        // 4. An SNV in the intron between that exon and the next (7,674,291-7,674,858).
        call(4, operation(0, 7674499, 7674501, "A", VariantType::SingleNucleotideVariant), None),
        // 5. An SNV past the last exon (both reads end at 7,687,490).
        call(5, operation(0, 7699999, 7700001, "A", VariantType::SingleNucleotideVariant), None),
        // 6. A 4-bp insertion two bases past the last exon end: within the insert length.
        call(6, operation(0, 7687492, 7687493, "ACGT", VariantType::Insertion), None),
        // 7. A 1-bp insertion at the same flank: one base out of reach.
        call(7, operation(0, 7687492, 7687493, "A", VariantType::Insertion), None),
        // 8. An SNV on chr18, where neither read has an exon.
        call(8, operation(1, 7674259, 7674261, "A", VariantType::SingleNucleotideVariant), None),
        // 9. A 1-bp insertion whose left flank is the first base of exon 7,674,181-7,674,290.
        call(9, operation(0, 7674181, 7674182, "A", VariantType::Insertion), None),
        // 10. The same one base to the left, on the last base of the intron before that exon.
        call(10, operation(0, 7674180, 7674181, "A", VariantType::Insertion), None)
    ];

    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = genotype_rna_variants(
        &read_ids,
        &variant_calls,
        &transcript_models_map
    );

    assert_eq!(genotypes.len(), variant_calls.len());
    let allele = |call_id: usize, read_id: usize| -> Allele {
        let alleles: &Vec<(usize, Allele)> = &genotypes[&call_id];
        assert_eq!(alleles.len(), read_ids.len());
        alleles.iter().find(|(id, _)| *id == read_id).unwrap().1.clone()
    };

    // 1. Member: Alternate. The reference read's exon spans the left flank: Reference.
    assert_eq!(allele(1, variant_read_id), Allele::Alternate);
    assert_eq!(allele(1, reference_read_id), Allele::Reference);
    // 2. The variant read's own insertion sits on this flank, so it is a sibling carrier:
    //    NotCovered, although its exon covers the site as the reference read's does.
    assert_eq!(allele(2, variant_read_id), Allele::NotCovered);
    assert_eq!(allele(2, reference_read_id), Allele::Reference);
    // 3. Both exons span the SNV.
    assert_eq!(allele(3, variant_read_id), Allele::Reference);
    assert_eq!(allele(3, reference_read_id), Allele::Reference);
    // 4. Spliced across.
    assert_eq!(allele(4, variant_read_id), Allele::NotCovered);
    assert_eq!(allele(4, reference_read_id), Allele::NotCovered);
    // 5. Ended before.
    assert_eq!(allele(5, variant_read_id), Allele::NotCovered);
    assert_eq!(allele(5, reference_read_id), Allele::NotCovered);
    // 6. A terminal clip may anchor up to its own length past the last aligned base.
    assert_eq!(allele(6, variant_read_id), Allele::Reference);
    assert_eq!(allele(6, reference_read_id), Allele::Reference);
    // 7. But not one base further.
    assert_eq!(allele(7, variant_read_id), Allele::NotCovered);
    assert_eq!(allele(7, reference_read_id), Allele::NotCovered);
    // 8. Another chromosome.
    assert_eq!(allele(8, variant_read_id), Allele::NotCovered);
    assert_eq!(allele(8, reference_read_id), Allele::NotCovered);
    // 9. An exon that starts on the flank holds it.
    assert_eq!(allele(9, variant_read_id), Allele::Reference);
    assert_eq!(allele(9, reference_read_id), Allele::Reference);
    // 10. An exon that starts one base past the flank does not, and the exon before it ended
    //     a whole intron earlier.
    assert_eq!(allele(10, variant_read_id), Allele::NotCovered);
    assert_eq!(allele(10, reference_read_id), Allele::NotCovered);
    // A read without a model is NotCovered everywhere.
    for call_id in 1..=10 {
        assert_eq!(allele(call_id, unmodelled_read_id), Allele::NotCovered);
    }
}

/// A breakend joining chr17:1,100 to chr17:900,000, and a translocation joining chr17:1,100 to
/// chr18:5,000. A read with an exon over 1,000 to 1,200 runs straight through side 1, and one with
/// an exon over 4,900 to 5,100 on chr18 through side 2 of the translocation: both hold the
/// reference there. No exon holds both sides of either call, so testing the span between them
/// counted no read as reference.
#[test]
fn site_span_tests_each_side_of_a_breakend_on_its_own() {
    let breakend = |chromosome_2: u16, position_2: u32, variant_type: VariantType| -> VariantCall {
        let mut variant_records: HashSet<VariantRecord> = HashSet::new();
        variant_records.insert(VariantRecord::new(
            1,
            100,
            101,
            GraphOperation::new(
                0,
                1_100,
                Strand::Forward,
                GraphOperationType::Downstream,
                chromosome_2,
                position_2,
                Strand::Forward,
                GraphOperationType::Upstream,
                "".into(),
                variant_type
            )
        ));
        VariantCall::from_variant_records(0, variant_records, 0, 4, 6, 2)
    };
    let breakpoint: SiteSpan = SiteSpan::from_variant_call(&breakend(0, 900_000, VariantType::Breakpoint));
    let translocation: SiteSpan = SiteSpan::from_variant_call(&breakend(1, 5_000, VariantType::Translocation));

    let through_side_1: Vec<(u16, u32, u32)> = vec![(0, 1_000, 1_200)];
    let through_side_2: Vec<(u16, u32, u32)> = vec![(1, 4_900, 5_100)];
    let elsewhere: Vec<(u16, u32, u32)> = vec![(0, 2_000, 3_000), (1, 6_000, 7_000)];

    assert!(breakpoint.is_covered_by(&through_side_1));
    assert!(!breakpoint.is_covered_by(&elsewhere));
    assert!(translocation.is_covered_by(&through_side_1));
    assert!(translocation.is_covered_by(&through_side_2));
    assert!(!translocation.is_covered_by(&elsewhere));
}
