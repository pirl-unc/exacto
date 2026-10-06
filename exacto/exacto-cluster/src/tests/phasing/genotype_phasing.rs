use std::collections::{HashMap, HashSet};
use exacto_caller::prelude::*;

use super::*;


/// A cell holds at least `min_reads` reads that tell it from every other cell.
///
/// 170 reads of one cluster are genotyped at two variant calls, a breakpoint (call 1) and an
/// insertion (call 2):
///
///   Reads     Breakpoint    Insertion
///   0-158     Alternate     Alternate
///   159-168   NotCovered    Alternate
///   169       Alternate     Reference
///
/// Read 169 carries the insertion spelled as two mismatches. Two cells of 169 reads and 1 read
/// leave no error and have the lowest BIC, 29.00 against 30.99 for one cell. Under a minimum
/// of 3 reads the two cells are not a choice, and the cluster is one cell. Under a minimum of
/// 1 read they are.
///
/// 170 reads of another cluster split 85 to 85 at the insertion, and the two cells stand under
/// the minimum of 3 reads.
///
/// 178 reads of a third cluster carry the breakpoint, and 3 of them carry the insertion spelled
/// as two mismatches. The cells of 175 reads and 3 reads stand under the minimum of 3 reads and
/// not under 21, the read support of a variant at a depth of 178 reads.
#[test]
fn phase_genotypes_returns_cells_with_min_reads() {
    let read_ids: HashSet<usize> = (0..170usize).collect();
    let mut breakpoint: Vec<(usize, Allele)> = Vec::new();
    let mut insertion: Vec<(usize, Allele)> = Vec::new();
    for read_id in 0..159usize {
        breakpoint.push((read_id, Allele::Alternate));
        insertion.push((read_id, Allele::Alternate));
    }
    for read_id in 159..169usize {
        breakpoint.push((read_id, Allele::NotCovered));
        insertion.push((read_id, Allele::Alternate));
    }
    breakpoint.push((169, Allele::Alternate));
    insertion.push((169, Allele::Reference));
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([(1, breakpoint), (2, insertion)]);

    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 3, 7, 1_000, 10_000, 42);
    let mut cell_read_ids: Vec<HashSet<usize>> = cells.into_values().collect();
    cell_read_ids.sort_by_key(|cell| std::cmp::Reverse(cell.len()));
    assert_eq!(cell_read_ids, vec![read_ids.clone()]);

    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 1, 7, 1_000, 10_000, 42);
    let mut cell_read_ids: Vec<HashSet<usize>> = cells.into_values().collect();
    cell_read_ids.sort_by_key(|cell| std::cmp::Reverse(cell.len()));
    assert_eq!(cell_read_ids, vec![(0..169usize).collect::<HashSet<usize>>(), HashSet::from([169])]);

    let insertion: Vec<(usize, Allele)> = (0..170usize)
        .map(|read_id| (read_id, if read_id < 85 { Allele::Alternate } else { Allele::Reference }))
        .collect();
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([(2, insertion)]);

    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 3, 7, 1_000, 10_000, 42);
    let mut cell_read_ids: Vec<HashSet<usize>> = cells.into_values().collect();
    cell_read_ids.sort_by_key(|cell| *cell.iter().min().unwrap());
    assert_eq!(
        cell_read_ids,
        vec![(0..85usize).collect::<HashSet<usize>>(), (85..170usize).collect::<HashSet<usize>>()]
    );

    let read_ids: HashSet<usize> = (0..178usize).collect();
    let breakpoint: Vec<(usize, Allele)> = (0..178usize).map(|read_id| (read_id, Allele::Alternate)).collect();
    let insertion: Vec<(usize, Allele)> = (0..178usize)
        .map(|read_id| (read_id, if read_id < 175 { Allele::Alternate } else { Allele::Reference }))
        .collect();
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([(1, breakpoint), (2, insertion)]);

    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 3, 7, 1_000, 10_000, 42);
    let mut cell_read_ids: Vec<HashSet<usize>> = cells.into_values().collect();
    cell_read_ids.sort_by_key(|cell| std::cmp::Reverse(cell.len()));
    assert_eq!(
        cell_read_ids,
        vec![(0..175usize).collect::<HashSet<usize>>(), (175..178usize).collect::<HashSet<usize>>()]
    );

    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 21, 7, 1_000, 10_000, 42);
    let cell_read_ids: Vec<HashSet<usize>> = cells.into_values().collect();
    assert_eq!(cell_read_ids, vec![read_ids.clone()]);

    // 297 reads of a fourth cluster are genotyped at two insertions:
    //
    //   Reads     First insertion   Second insertion
    //   0-104     Alternate         Reference
    //   105-251   Reference         Reference
    //   252-259   Reference         Alternate
    //   260-286   Reference         NotCovered
    //   287-296   NotCovered        NotCovered
    //
    // Reads 260-286 fit the cell of reads 105-251 and the cell of reads 252-259 alike, and
    // reads 287-296 fit every cell.
    let read_ids: HashSet<usize> = (0..297usize).collect();
    let alleles = |read_id: usize| -> (Allele, Allele) {
        match read_id {
            0..=104 => (Allele::Alternate, Allele::Reference),
            105..=251 => (Allele::Reference, Allele::Reference),
            252..=259 => (Allele::Reference, Allele::Alternate),
            260..=286 => (Allele::Reference, Allele::NotCovered),
            _ => (Allele::NotCovered, Allele::NotCovered)
        }
    };
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([
        (1, (0..297usize).map(|read_id| (read_id, alleles(read_id).0)).collect()),
        (2, (0..297usize).map(|read_id| (read_id, alleles(read_id).1)).collect())
    ]);

    // Under a minimum of 1 read and this seed, both are assigned to the cell of reads 252-259, a
    // cell of 45 reads. Each is shared with the other cells it fits: reads 260-286 with the cell
    // of reads 105-251, reads 287-296 with both other cells.
    let (cells, shared) = phase_genotypes(&read_ids, &genotypes, 1, 7, 1_000, 10_000, 9);
    let mut num_reads: Vec<usize> = cells.values().map(|cell| cell.len()).collect();
    num_reads.sort_unstable();
    assert_eq!(num_reads, vec![45, 105, 147]);
    assert_eq!(cells.values().find(|cell| cell.contains(&252)).unwrap().len(), 45);
    let shared_of = |read_id: usize| -> HashSet<usize> {
        let cell: usize = *cells.iter().find(|(_, reads)| reads.contains(&read_id)).unwrap().0;
        shared.get(&cell).cloned().unwrap_or_default()
    };
    assert_eq!(shared_of(0), (287..297usize).collect::<HashSet<usize>>());
    assert_eq!(shared_of(105), (260..297usize).collect::<HashSet<usize>>());
    assert_eq!(shared_of(252), HashSet::new());

    // The reads telling that cell from the others are 8, under a minimum of 31. The cluster is
    // two cells, split at the first insertion, and reads 287-296 are assigned to either and
    // shared with the other.
    let (cells, shared) = phase_genotypes(&read_ids, &genotypes, 31, 7, 1_000, 10_000, 9);
    assert_eq!(cells.len(), 2);
    let carriers: &HashSet<usize> = cells.values().find(|cell| cell.contains(&0)).unwrap();
    let others: &HashSet<usize> = cells.values().find(|cell| cell.contains(&105)).unwrap();
    assert!((0..105usize).all(|read_id| carriers.contains(&read_id)));
    assert!((105..287usize).all(|read_id| others.contains(&read_id)));
    for (cell, reads) in cells.iter() {
        let unassigned: HashSet<usize> = (287..297usize).filter(|read_id| !reads.contains(read_id)).collect();
        assert_eq!(shared.get(cell).cloned().unwrap_or_default(), unassigned);
    }
}


/// A read is also a read of every other cell it does not tell its own cell from, and of no
/// other.
///
///   Reads     Site 1        Site 2
///   0-39      Alternate     Alternate
///   40-79     Reference     Reference
///   80        Alternate     Reference
///   81-85     NotCovered    NotCovered
///
/// Under a minimum of 3 reads the cluster is two cells, of reads 0-39 and of reads 40-79. Reads
/// 81-85 cover neither site, so each fits both cells: it is assigned to one and shared with the
/// other. Read 80 fits each cell with one error, but it tells the two apart at both sites, so it
/// is in one cell only: a read at odds with both cells is evidence of an error, not of both.
/// Under a minimum of 41 reads no site is phased, the cluster is one cell, and nothing is
/// shared.
#[test]
fn phase_genotypes_shares_reads_that_do_not_tell_cells_apart() {
    let read_ids: HashSet<usize> = (0..86usize).collect();
    let alleles = |read: usize| -> (Allele, Allele) {
        match read {
            0..=39 => (Allele::Alternate, Allele::Alternate),
            40..=79 => (Allele::Reference, Allele::Reference),
            80 => (Allele::Alternate, Allele::Reference),
            _ => (Allele::NotCovered, Allele::NotCovered)
        }
    };
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([
        (1, (0..86usize).map(|read| (read, alleles(read).0)).collect()),
        (2, (0..86usize).map(|read| (read, alleles(read).1)).collect())
    ]);

    let (cells, shared) = phase_genotypes(&read_ids, &genotypes, 3, 7, 1_000, 10_000, 42);
    assert_eq!(cells.len(), 2);
    let cell_of = |read: usize| -> usize { *cells.iter().find(|(_, reads)| reads.contains(&read)).unwrap().0 };
    assert!((0..40).all(|read| cell_of(read) == cell_of(0)));
    assert!((40..80).all(|read| cell_of(read) == cell_of(40)));
    assert_ne!(cell_of(0), cell_of(40));
    for (cell, reads) in cells.iter() {
        let unassigned: HashSet<usize> = (81..86usize).filter(|read| !reads.contains(read)).collect();
        assert_eq!(shared.get(cell).cloned().unwrap_or_default(), unassigned, "cell {cell}");
    }
    assert!(shared.values().all(|reads| !reads.contains(&80)));

    let (cells, shared) = phase_genotypes(&read_ids, &genotypes, 41, 7, 1_000, 10_000, 42);
    assert_eq!(cells.into_values().collect::<Vec<HashSet<usize>>>(), vec![read_ids]);
    assert!(shared.is_empty());
}


#[test]
fn phase_genotypes_keeps_sites_short_of_the_floor_from_moving_reads() {
    for relabelling in 0..8usize {
        // A read's id under this relabelling: a bijection of 0..220.
        let id = |read: usize| -> usize { (read * 97 + relabelling * 31) % 220 };
        let read_ids: HashSet<usize> = (0..220usize).map(id).collect();
        let mut site_1: Vec<(usize, Allele)> = Vec::new();
        let mut site_2: Vec<(usize, Allele)> = Vec::new();
        let mut site_3: Vec<(usize, Allele)> = Vec::new();
        for read in 0..220usize {
            let (allele_1, allele_2): (Allele, Allele) = match read {
                0..=99 => (Allele::Alternate, if read < 30 { Allele::Reference } else { Allele::NotCovered }),
                100..=199 => (Allele::Reference, if read < 130 { Allele::Reference } else { Allele::NotCovered }),
                _ => (Allele::NotCovered, Allele::Alternate)
            };
            site_1.push((id(read), allele_1));
            site_2.push((id(read), allele_2.clone()));
            site_3.push((id(read), allele_2));
        }
        let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([(1, site_1), (2, site_2), (3, site_3)]);
        let (cells, _) = phase_genotypes(&read_ids, &genotypes, 35, 7, 1_000, 10_000, 42);
        assert_eq!(cells.len(), 2, "relabelling {relabelling}");
        let cell_of = |read: usize| -> usize { *cells.iter().find(|(_, reads)| reads.contains(&id(read))).unwrap().0 };
        assert!((0..100).all(|read| cell_of(read) == cell_of(0)), "relabelling {relabelling}: A split");
        assert!((100..200).all(|read| cell_of(read) == cell_of(100)), "relabelling {relabelling}: B split");
        assert!((200..220).all(|read| cell_of(read) == cell_of(200)), "relabelling {relabelling}: C split");
        assert_ne!(cell_of(0), cell_of(100), "relabelling {relabelling}");
    }
}


#[test]
fn phase_genotypes_returns_one_cell_for_sites_short_of_the_floor() {
    let read_ids: HashSet<usize> = (0..285usize).collect();
    let alleles = |read: usize| -> (Allele, Allele) {
        match read {
            0..=199 => (Allele::Alternate, Allele::NotCovered),
            200..=204 => (Allele::Reference, Allele::Reference),
            205..=264 => (Allele::NotCovered, Allele::Reference),
            _ => (Allele::NotCovered, Allele::Alternate)
        }
    };
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([
        (1, (0..285usize).map(|read| (read, alleles(read).0)).collect()),
        (2, (0..285usize).map(|read| (read, alleles(read).1)).collect())
    ]);
    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 44, 7, 1_000, 10_000, 42);
    assert_eq!(cells.into_values().collect::<Vec<HashSet<usize>>>(), vec![read_ids]);
}


#[test]
fn phase_genotypes_returns_one_cell_below_twice_the_floor() {
    let read_ids: HashSet<usize> = (0..60usize).collect();
    let site: Vec<(usize, Allele)> = (0..60usize)
        .map(|read| (read, if read < 30 { Allele::Alternate } else { Allele::Reference }))
        .collect();
    let genotypes: HashMap<usize, Vec<(usize, Allele)>> = HashMap::from([(1, site)]);
    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 31, 7, 1_000, 10_000, 42);
    assert_eq!(cells.into_values().collect::<Vec<HashSet<usize>>>(), vec![read_ids.clone()]);
    let (cells, _) = phase_genotypes(&read_ids, &genotypes, 30, 7, 1_000, 10_000, 42);
    assert_eq!(cells.len(), 2);
}
