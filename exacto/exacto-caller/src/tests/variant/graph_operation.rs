use super::*;


/// Cryptic exons, retained introns and UTR extensions are one block of bases, typed after it is
/// built (`calling/rna/variant_record_caller.rs`), and both of its positions are bases of the block.
/// rna-008's cryptic exon, chr17:7672205-7672484, is 280 bases.
#[test]
fn get_variant_size_counts_both_ends_of_a_block() {
    for variant_type in [VariantType::CrypticExon, VariantType::IntronRetention, VariantType::UTRExtension] {
        for (position_2, size) in [(7_672_484, 280), (7_672_205, 1)] {
            let block: GraphOperation = GraphOperation::new(
                0,
                7_672_205,
                Strand::Reverse,
                GraphOperationType::Include,
                0,
                position_2,
                Strand::Reverse,
                GraphOperationType::Include,
                "".into(),
                variant_type.clone()
            );
            assert_eq!(block.get_variant_size(), size, "{} ending at {}", variant_type.as_str(), position_2);
        }
    }
}
