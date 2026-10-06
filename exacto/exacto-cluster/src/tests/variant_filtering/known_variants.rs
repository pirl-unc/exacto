use exacto_core::prelude::Strand;

use super::*;


/// A call matches an allowed variant with the same type and chromosomes, both positions within
/// the distance, and an insertion's sequence within the edit distance.
///
/// Listed: a 12-base insertion at chr17:7,674,224-7,674,225 and a breakend 7,669,660-7,673,701.
///
///   Call                                      Max distance   Allowed
///   the insertion as listed                   10             yes
///   the insertion 10 bases on                 10             yes
///   the insertion 11 bases on                 10             no
///   the insertion with a base more            10             yes, 1 edit in 13 bases
///   another sequence at the same place        10             no
///   a deletion at the insertion's positions   10             no, another type
///   the insertion on chr18                    10             no
///   the breakend 9 bases off                  10             yes
///   the breakend 9 bases off                  5              no
#[test]
fn allows_matches_type_chromosomes_positions_and_sequence() {
    let known_variants: KnownVariants = KnownVariants {
        dna_variants: HashSet::new(),
        allowed_variants: Some(vec![
            (0, 0, 7_674_224, 7_674_225, VariantType::Insertion, "CCCATCCGCCTG".into()),
            (0, 0, 7_669_660, 7_673_701, VariantType::Breakpoint, "".into())
        ])
    };
    let op = |chromosome: u16, position_1: u32, position_2: u32, sequence: &str, variant_type: VariantType| -> GraphOperation {
        GraphOperation::new(
            chromosome, position_1, Strand::Reverse, GraphOperationType::Downstream,
            chromosome, position_2, Strand::Reverse, GraphOperationType::Upstream,
            sequence.into(), variant_type
        )
    };
    // (call, max distance, allowed)
    let cases: Vec<(GraphOperation, u32, bool)> = vec![
        (op(0, 7_674_224, 7_674_225, "CCCATCCGCCTG", VariantType::Insertion), 10, true),
        (op(0, 7_674_234, 7_674_235, "CCCATCCGCCTG", VariantType::Insertion), 10, true),
        (op(0, 7_674_235, 7_674_236, "CCCATCCGCCTG", VariantType::Insertion), 10, false),
        (op(0, 7_674_224, 7_674_225, "CCCATCCGCCCTG", VariantType::Insertion), 10, true),
        (op(0, 7_674_224, 7_674_225, "AATTAATTAATT", VariantType::Insertion), 10, false),
        (op(0, 7_674_224, 7_674_225, "", VariantType::Deletion), 10, false),
        (op(1, 7_674_224, 7_674_225, "CCCATCCGCCTG", VariantType::Insertion), 10, false),
        (op(0, 7_669_651, 7_673_701, "", VariantType::Breakpoint), 10, true),
        (op(0, 7_669_651, 7_673_701, "", VariantType::Breakpoint), 5, false)
    ];
    for (call, max_distance, expected) in cases {
        assert_eq!(
            known_variants.allows(&call, max_distance, 0.5),
            expected,
            "{:?} {}-{} {}, max distance {}",
            call.get_variant_type(),
            call.get_position_1(),
            call.get_position_2(),
            call.get_sequence(),
            max_distance
        );
    }
}
