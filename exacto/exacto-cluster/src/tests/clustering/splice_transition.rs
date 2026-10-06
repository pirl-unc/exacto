use exacto_caller::prelude::{ReadFilter, RNAJunctionReadSupportFilter, RNAJunctionReadSupportIndex};
use exacto_core::prelude::Strand;

use super::*;


/// An unannotated transition is held to the shared splicing event error model, even when every
/// junction of the chain is supported on its own.
///
/// Three junctions a, b and c of one chromosome follow one another in a reference transcript
/// (a -> b and b -> c are annotated):
///
///   286 reads  a b c a b c
///    14 reads  a c a b c
///
/// a -> c is spliced by the 14 reads of the 300 that leave a or enter c, against a minimum of 46
/// at an error rate of 0.01 and a max FPR of 1e-6, so their chain is dropped. c -> a is spliced
/// by all 300.
#[test]
fn splice_transition_read_counts_fails_chain_with_unsupported_order_of_supported_junctions() {
    let a: SpliceJunction = SpliceJunction::new(0, 0, 100, 150, Strand::Forward, Strand::Forward);
    let b: SpliceJunction = SpliceJunction::new(0, 0, 300, 350, Strand::Forward, Strand::Forward);
    let c: SpliceJunction = SpliceJunction::new(0, 0, 500, 550, Strand::Forward, Strand::Forward);
    let good: Vec<SpliceJunction> = vec![a.clone(), b.clone(), c.clone(), a.clone(), b.clone(), c.clone()];
    let bad: Vec<SpliceJunction> = vec![a.clone(), c.clone(), a.clone(), b.clone(), c.clone()];
    let chains: HashMap<Vec<SpliceJunction>, HashSet<usize>> = HashMap::from([
        (good.clone(), (0..286).collect()),
        (bad.clone(), (286..300).collect())
    ]);
    let annotated_transitions: HashSet<(SpliceJunction, SpliceJunction)> = HashSet::from([
        (a.clone(), b.clone()),
        (b.clone(), c.clone())
    ]);

    let read_counts: SpliceTransitionReadCounts = SpliceTransitionReadCounts::new(&chains, &annotated_transitions, &HashMap::new());
    let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(&read_counts.get_depths(), 0.01, 1e-6, 1);
    let read_filter: RNAJunctionReadSupportFilter<SpliceTransitionReadCounts> = RNAJunctionReadSupportFilter::new(&read_counts, &read_support_index);
    let failed_chains: HashSet<Vec<SpliceJunction>> = chains.keys().filter(|chain| !read_filter.passes(chain)).cloned().collect();

    assert_eq!(read_counts.get_read_counts(&good), vec![(300, 300)]);
    assert_eq!(read_counts.get_read_counts(&bad), vec![(14, 300), (300, 300)]);
    assert_eq!(failed_chains, HashSet::from([bad]));
}


/// A supported novel transition and an annotated one pass, and the laps of one read count it
/// once.
///
/// Junctions a, b and c of one chromosome, with a -> b annotated:
///
///   250 reads  a b
///    50 reads  a c
///
/// a -> c is spliced by 50 of 300 reads and passes. In place of the 50 reads, 14 reads spell
/// a c a c a c: a -> c is spliced by 14 of the 264 that leave a or enter c, short of its minimum,
/// and their chain fails. With every transition annotated nothing is tested.
#[test]
fn splice_transition_read_counts_passes_supported_and_annotated_transitions_counting_laps_once() {
    let a: SpliceJunction = SpliceJunction::new(0, 0, 100, 150, Strand::Forward, Strand::Forward);
    let b: SpliceJunction = SpliceJunction::new(0, 0, 300, 350, Strand::Forward, Strand::Forward);
    let c: SpliceJunction = SpliceJunction::new(0, 0, 500, 550, Strand::Forward, Strand::Forward);
    let ordinary: Vec<SpliceJunction> = vec![a.clone(), b.clone()];
    let novel: Vec<SpliceJunction> = vec![a.clone(), c.clone()];
    let repeats: Vec<SpliceJunction> = vec![a.clone(), c.clone(), a.clone(), c.clone(), a.clone(), c.clone()];
    let annotated_transitions: HashSet<(SpliceJunction, SpliceJunction)> = HashSet::from([(a.clone(), b.clone())]);
    let all_annotated_transitions: HashSet<(SpliceJunction, SpliceJunction)> = HashSet::from([
        (a.clone(), b.clone()),
        (a.clone(), c.clone()),
        (c.clone(), a.clone())
    ]);

    // (chains, annotated transitions, chains that fail)
    let cases: Vec<(HashMap<Vec<SpliceJunction>, HashSet<usize>>, &HashSet<(SpliceJunction, SpliceJunction)>, HashSet<Vec<SpliceJunction>>)> = vec![
        (
            HashMap::from([(ordinary.clone(), (0..250).collect()), (novel.clone(), (250..300).collect())]),
            &annotated_transitions,
            HashSet::new()
        ),
        (
            HashMap::from([(ordinary.clone(), (0..250).collect()), (repeats.clone(), (250..264).collect())]),
            &annotated_transitions,
            HashSet::from([repeats.clone()])
        ),
        (
            HashMap::from([(ordinary.clone(), (0..250).collect()), (repeats.clone(), (250..264).collect())]),
            &all_annotated_transitions,
            HashSet::new()
        )
    ];
    for (chains, annotated_transitions, expected_failed_chains) in cases {
        let read_counts: SpliceTransitionReadCounts = SpliceTransitionReadCounts::new(&chains, annotated_transitions, &HashMap::new());
        let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(&read_counts.get_depths(), 0.01, 1e-6, 1);
        let read_filter: RNAJunctionReadSupportFilter<SpliceTransitionReadCounts> = RNAJunctionReadSupportFilter::new(&read_counts, &read_support_index);
        let failed_chains: HashSet<Vec<SpliceJunction>> = chains.keys().filter(|chain| !read_filter.passes(chain)).cloned().collect();

        assert_eq!(failed_chains, expected_failed_chains);
    }
}
