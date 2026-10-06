use exacto_core::prelude::*;
use noodles_bam as bam;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use super::*;


/// Editing keeps the model whole: `remove` and `retain` rebase every read position that bases,
/// events and records carry, and the quality cache follows the bases.
///
/// `scga-mini-rna-013` read 182 has a Reverse primary (`6S`, two introns, seven deletions and a
/// `1I`, then `459S`) and a Forward supplementary with five introns, `1180I/13I/1I` and two `1D`,
/// so events and records sit on both strands.
#[test]
fn editing_remove_and_retain_keep_the_model_whole() {
    let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_file.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-013-tumor_chunk_0000/182/ccs").unwrap();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
    assert_eq!(original.num_bases(), 1900);
    assert_eq!(original.get_records().len(), 2);
    assert!(original.num_events() > 0);

    // An event is the same event when its ends stand on the same placed bases.
    let signature = |model: &AlignmentModel, (previous, next): (u32, u32)| -> String {
        format!(
            "{:?} {:?} {}",
            model.get_base(previous).get_placement(),
            model.get_base(next).get_placement(),
            model.get_event(previous, next).unwrap().get_kind().as_str()
        )
    };

    // Cut ten bases out of the middle. Whatever they carried, the rest must line up.
    let range = 1000..1010u32;
    let mut model: AlignmentModel = original.clone();
    let mut editor = AlignmentModelEditor::new(&mut model);
    editor.remove(range.clone());
    assert_eq!(editor.get_history(), &[AlignmentModelEdit::Removal { range: range.clone() }]);
    assert!((0..editor.get_model().num_bases()).all(|position| !editor.is_restored(position)), "removing restores nothing");
    assert_eq!(model.get_read_sequence(), format!("{}{}", &sequence[..1000], &sequence[1010..]));
    assert_eq!(model.iter_base_quality_scores().collect::<Vec<_>>(), &[&qualities[..1000], &qualities[1010..]].concat()[..]);
    assert_eq!(model.iter_base_quality_scores().collect::<Vec<_>>(), model.get_bases().iter().map(AlignmentModelBase::get_base_quality).collect::<Vec<u8>>());
    for (position, base) in model.get_bases().iter().enumerate() {
        assert_eq!(base.get_read_position(), position as u32);
    }
    let mut before: Vec<String> = original.get_events().keys()
        .filter(|(previous, next)| !range.contains(previous) && !range.contains(next))
        .map(|&key| signature(&original, key))
        .collect();
    let mut after: Vec<String> = model.get_events().keys().map(|&key| signature(&model, key)).collect();
    before.sort();
    after.sort();
    assert_eq!(after, before, "the surviving events are the same events, found at their bases' new positions");
    assert!(model.num_events() < original.num_events() || before.len() == original.get_events().len());
    for (&(previous, next), event) in model.get_events() {
        assert_eq!((event.get_prev_read_position(), event.get_next_read_position()), (previous, next));
        assert!(next < model.num_bases());
        assert!(model.get_event_at(previous).is_some() && model.get_event_at(next).is_some(), "the index mirrors the keys");
    }
    assert_eq!(model.get_records().len(), 2);
    for (record, before) in model.get_records().iter().zip(original.get_records()) {
        assert!(record.read_start <= record.read_end && record.read_end < model.num_bases());
        assert_eq!(model.get_base(record.read_start).get_placement(), original.get_base(before.read_start).get_placement());
        assert_eq!(model.get_base(record.read_end).get_placement(), original.get_base(before.read_end).get_placement());
        let spanned: u32 = range.end.min(before.read_end + 1).saturating_sub(range.start.max(before.read_start));
        assert_eq!(
            (before.read_end - before.read_start) - (record.read_end - record.read_start),
            spanned,
            "a record shrinks by exactly the part of the cut it spanned"
        );
    }

    // Keep a middle range that ends on a parked deletion: the cut takes the deletion with it.
    let cut: u32 = original.get_bases().iter().filter_map(AlignmentModelBase::get_deletion_read_position).max().unwrap();
    let mut model: AlignmentModel = original.clone();
    let mut editor = AlignmentModelEditor::new(&mut model);
    editor.retain(500..cut);
    assert_eq!(editor.get_history(), &[AlignmentModelEdit::Retention { range: 500..cut }]);
    assert!((0..editor.get_model().num_bases()).all(|position| !editor.is_restored(position)), "retaining restores nothing");
    assert_eq!(model.num_bases(), cut - 500);
    assert_eq!(model.get_read_sequence(), &sequence[500..cut as usize]);
    assert_eq!(model.iter_base_quality_scores().collect::<Vec<_>>(), &qualities[500..cut as usize]);
    for (position, base) in model.get_bases().iter().enumerate() {
        assert_eq!(base.get_read_position(), position as u32);
        assert!(base.get_deletion_read_position().is_none_or(|position| position < model.num_bases()), "no deletion is parked at the cut");
    }
    let mut before: Vec<String> = original.get_events().keys()
        .filter(|(previous, next)| (500..cut).contains(previous) && (500..cut).contains(next))
        .map(|&key| signature(&original, key))
        .collect();
    let mut after: Vec<String> = model.get_events().keys().map(|&key| signature(&model, key)).collect();
    before.sort();
    after.sort();
    assert_eq!(after, before);
    for (&(previous, next), _) in model.get_events() {
        assert!(model.get_event_at(previous).is_some() && model.get_event_at(next).is_some());
    }
    for record in model.get_records() {
        assert!(record.read_start <= record.read_end && record.read_end < model.num_bases());
    }
}


/// The base-level edits: a parked deletion comes back on its flank, placed and in read
/// orientation on either strand; a substitution follows the reference; unplaced bases are soft
/// clips. Each leaves the quality cache matching the bases.
///
/// `scga-mini-rna-007` read 569 (Forward, `195=1I570=1D439=1I1730=1I725=`) and `scga-mini-rna-001`
/// read 33 (Reverse, `416=1D274=`) are single-record reads that each park exactly one deletion,
/// cross no intron and carry no mismatch: the first such read in each file. No scga-mini file
/// holds one on each strand.
#[test]
fn editing_restores_substitutes_and_inserts_on_the_models_own_bases() {
    for (bam_name, read_name, expected_strand) in [
        ("alignment/scga-mini-rna-007-tumor_minimap2_sorted.bam", "scga-mini-rna-007-tumor_chunk_0000/569/ccs", Strand::Forward),
        ("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam", "scga-mini-rna-001-tumor_chunk_0000/33/ccs", Strand::Reverse)
    ] {
        let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(bam_name)).unwrap();
        let bam_file: &str = bam_file.to_str().unwrap();
        let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
        let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
        reader.read_header().unwrap();
        let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
        let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
        let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
        let mut model: AlignmentModel = original.clone();

        // The one parked deletion: its flank, its payload and where the flank stands.
        let flank: &AlignmentModelBase = original.get_bases().iter().find(|base| base.get_deletion_read_position().is_some()).unwrap();
        assert_eq!(original.get_bases().iter().filter(|base| base.get_deletion_read_position().is_some()).count(), 1);
        assert_eq!(original.num_events(), 1);
        let position: u32 = flank.get_deletion_read_position().unwrap();
        let payload: Vec<Nucleotide> = flank.get_deleted_reference_bases().to_vec();
        let (chromosome_id, reference, strand) = flank.get_placement().get_coordinate().unwrap();
        let strand: Strand = strand.clone();
        assert_eq!(strand, expected_strand);
        let flank: u32 = flank.get_read_position();
        assert_eq!(position, if strand == Strand::Forward { flank + 1 } else { flank });
        let length: u32 = payload.len() as u32;

        let deletion: u32 = position;
        let mut editor = AlignmentModelEditor::new(&mut model);
        assert!(editor.restore_deletion(position, 7));
        assert_eq!(editor.get_model().num_bases(), original.num_bases() + length);
        for (index, nucleotide) in payload.iter().enumerate() {
            let base: &AlignmentModelBase = editor.get_model().get_base(position + index as u32);
            assert!(editor.is_restored(position + index as u32));
            assert_eq!(base.get_nucleotide(), nucleotide);
            assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Match);
            assert_eq!(base.get_base_quality(), 7);
            let expected: u32 = if strand == Strand::Forward { reference + index as u32 + 1 } else { reference + length - index as u32 };
            assert_eq!(base.get_placement().get_coordinate(), Some((chromosome_id, expected, &strand)), "restored bases run away from the flank along the strand");
        }
        assert_eq!((0..editor.get_model().num_bases()).filter(|&position| editor.is_restored(position)).count() as u32, length, "only the put-back bases are restored");
        // The flank stays put on the forward strand and moves past the restored bases on the
        // reverse strand; either way it parks nothing now and its event is gone.
        let flank_now: u32 = if strand == Strand::Forward { flank } else { flank + length };
        assert_eq!(editor.get_model().get_base(flank_now).get_placement().get_coordinate(), Some((chromosome_id, reference, &strand)));
        assert!(editor.get_model().get_base(flank_now).get_deleted_reference_bases().is_empty());
        assert_eq!(editor.get_model().num_events(), 0);
        assert!(!editor.restore_deletion(position, 7), "nothing is parked there any more");
        assert_eq!(editor.get_model().get_read_sequence().len() as u32, editor.get_model().num_bases());
        assert_eq!(editor.get_model().iter_base_quality_scores().collect::<Vec<_>>(), editor.get_model().get_bases().iter().map(AlignmentModelBase::get_base_quality).collect::<Vec<u8>>());
        for (index, base) in editor.get_model().get_bases().iter().enumerate() {
            assert_eq!(base.get_read_position(), index as u32);
        }

        // A substitution follows the reference: rewritten away from it the base is a mismatch,
        // and back again a match. The same base is no change.
        let position: u32 = (0..editor.get_model().num_bases())
            .find(|&position| *editor.get_model().get_base(position).get_kind() == AlignmentModelBaseKind::Match && !editor.is_restored(position))
            .unwrap();
        let sequenced: Nucleotide = editor.get_model().get_base(position).get_nucleotide().clone();
        let other: Nucleotide = if sequenced == Nucleotide::A { Nucleotide::C } else { Nucleotide::A };
        assert!(!editor.substitute(position, sequenced.clone(), 9));
        assert!(editor.substitute(position, other.clone(), 9));
        assert_eq!(*editor.get_model().get_base(position).get_kind(), AlignmentModelBaseKind::Mismatch);
        assert_eq!(editor.get_model().get_base(position).get_reference_nucleotide(), Some(&sequenced));
        assert_eq!(editor.get_model().iter_base_quality_scores().collect::<Vec<_>>()[position as usize], 9);
        assert!(editor.substitute(position, sequenced.clone(), 11));
        assert_eq!(*editor.get_model().get_base(position).get_kind(), AlignmentModelBaseKind::Match);
        assert_eq!(editor.get_model().iter_base_quality_scores().collect::<Vec<_>>()[position as usize], 11);
        assert_eq!(editor.get_model().num_bases(), original.num_bases() + length);
        assert!(!editor.is_restored(position), "a substitution rewrites a sequenced base; it adds none");

        // Unplaced bases are soft clips that carry no placement, in front of the base they name.
        let before: u32 = editor.get_model().num_bases();
        editor.insert_unplaced(position, &[Nucleotide::G, Nucleotide::T], 3);
        assert_eq!(editor.get_model().num_bases(), before + 2);
        for offset in 0..2 {
            let base: &AlignmentModelBase = editor.get_model().get_base(position + offset);
            assert!(editor.is_restored(position + offset));
            assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Softclip);
            assert!(!base.get_placement().is_placed());
            assert_eq!(base.get_base_quality(), 3);
        }
        assert!(!editor.is_restored(position + 2));
        assert_eq!(&editor.get_model().get_read_sequence()[position as usize..position as usize + 2], "GT");
        assert_eq!(editor.get_model().get_base(position + 2).get_nucleotide(), &sequenced);
        assert_eq!(&editor.get_model().iter_base_quality_scores().collect::<Vec<_>>()[position as usize..position as usize + 3], &[3, 3, 11]);
        for (index, base) in editor.get_model().get_bases().iter().enumerate() {
            assert_eq!(base.get_read_position(), index as u32);
        }

        // The history holds the four edits that changed the model, oldest first; the refused
        // restoration and the no-change substitution are not there.
        assert_eq!(editor.get_history(), &[
            AlignmentModelEdit::DeletionRestoration { read_position: deletion, quality: 7 },
            AlignmentModelEdit::Substitution { read_position: position, nucleotide: other.clone(), quality: 9 },
            AlignmentModelEdit::Substitution { read_position: position, nucleotide: sequenced.clone(), quality: 11 },
            AlignmentModelEdit::UnplacedInsertion { read_position: position, nucleotides: vec![Nucleotide::G, Nucleotide::T], quality: 3 }
        ]);

        // Replaying the history on the starting model reproduces the edited model and its
        // restored bases.
        let history: Vec<AlignmentModelEdit> = editor.get_history().to_vec();
        let restored: Vec<bool> = (0..editor.get_model().num_bases()).map(|position| editor.is_restored(position)).collect();
        let mut replayed: AlignmentModel = original.clone();
        let mut replay = AlignmentModelEditor::new(&mut replayed);
        for edit in history {
            match edit {
                AlignmentModelEdit::Substitution { read_position, nucleotide, quality } => assert!(replay.substitute(read_position, nucleotide, quality)),
                AlignmentModelEdit::Removal { range } => replay.remove(range),
                AlignmentModelEdit::DeletionRestoration { read_position, quality } => assert!(replay.restore_deletion(read_position, quality)),
                AlignmentModelEdit::UnplacedInsertion { read_position, nucleotides, quality } => replay.insert_unplaced(read_position, &nucleotides, quality),
                AlignmentModelEdit::Retention { range } => replay.retain(range)
            }
        }
        assert_eq!((0..replay.get_model().num_bases()).map(|position| replay.is_restored(position)).collect::<Vec<bool>>(), restored);
        assert_eq!(replay.get_history(), editor.get_history());
        assert_eq!(format!("{:?}", replay.get_model().get_bases()), format!("{:?}", editor.get_model().get_bases()));
        assert_eq!(replay.get_model().get_events().len(), editor.get_model().get_events().len());
        for (record, edited) in replay.get_model().get_records().iter().zip(editor.get_model().get_records()) {
            assert_eq!((record.read_start, record.read_end), (edited.read_start, edited.read_end));
        }
    }
}


/// A substitution follows the reference only on aligned bases. An inserted base and a clipped
/// base keep their kind and placement and gain no reference nucleotide; only the nucleotide and
/// its quality change.
///
/// `scga-mini-rna-013` read 182: 282..=458 are the head of the supplementary's `1180I`, placed at
/// its anchor, and 1894..=1899 are the primary's `6S`.
#[test]
fn substitute_keeps_the_kind_and_placement_of_bases_off_the_reference() {
    let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_file.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-013-tumor_chunk_0000/182/ccs").unwrap();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
    let mut model: AlignmentModel = original.clone();

    for (position, kind) in [(300u32, AlignmentModelBaseKind::Insertion), (1899, AlignmentModelBaseKind::Softclip)] {
        let base: &AlignmentModelBase = original.get_base(position);
        assert_eq!(*base.get_kind(), kind);
        assert!(base.get_placement().is_placed());
        assert_eq!(base.get_reference_nucleotide(), None);
        let other: Nucleotide = if *base.get_nucleotide() == Nucleotide::A { Nucleotide::C } else { Nucleotide::A };

        assert!(AlignmentModelEditor::new(&mut model).substitute(position, other.clone(), 5));
        let edited: &AlignmentModelBase = model.get_base(position);
        assert_eq!(edited.get_nucleotide(), &other);
        assert_eq!(edited.get_base_quality(), 5);
        assert_eq!(*edited.get_kind(), kind, "the base keeps its kind");
        assert_eq!(edited.get_placement(), base.get_placement(), "and its placement");
        assert_eq!(edited.get_reference_nucleotide(), None, "and gains no reference nucleotide");
    }
    assert_eq!(model.num_bases(), original.num_bases());
    assert_eq!(model.num_events(), original.num_events());
}


/// A cut that takes every base of a record drops the record, whether the cut starts at the read
/// start or further in. The other record is clipped and rebased, only the events with both ends
/// outside the cut survive, and parked deletions travel with their flanks.
///
/// `scga-mini-rna-013` read 182: the Forward supplementary covers 0..=1563 and the Reverse primary
/// 459..=1893. The primary parks its last two deletions at 1657 and 1851.
#[test]
fn remove_drops_a_record_with_no_base_left() {
    let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_file.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-013-tumor_chunk_0000/182/ccs").unwrap();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
    let spans = |model: &AlignmentModel| -> Vec<(u32, u32, Strand)> {
        model.get_records().iter().map(|record| (record.read_start, record.read_end, record.reference_strand.clone())).collect()
    };
    let events = |model: &AlignmentModel| -> Vec<((u32, u32), AlignmentModelEventKind)> {
        let mut events: Vec<((u32, u32), AlignmentModelEventKind)> = model.get_events().iter().map(|(&key, event)| (key, event.get_kind().clone())).collect();
        events.sort_by_key(|(key, _)| *key);
        events
    };
    assert_eq!(spans(&original), [(0, 1563, Strand::Forward), (459, 1893, Strand::Reverse)]);

    // From the read start through the supplementary's last base: the primary's tail is left.
    let mut model: AlignmentModel = original.clone();
    AlignmentModelEditor::new(&mut model).remove(0..1564);
    assert_eq!(model.get_read_sequence(), &sequence[1564..]);
    assert_eq!(spans(&model), [(0, 329, Strand::Reverse)]);
    assert_eq!(events(&model), [((92, 93), AlignmentModelEventKind::Deletion), ((286, 287), AlignmentModelEventKind::Deletion)]);
    assert_eq!(model.get_bases().iter().filter_map(AlignmentModelBase::get_deletion_read_position).collect::<Vec<u32>>(), [93, 287]);
    for &(previous, next) in model.get_events().keys() {
        assert_eq!(model.get_event_at(previous).map(AlignmentModelEvent::get_next_read_position), Some(next));
        assert_eq!(model.get_event_at(next).map(AlignmentModelEvent::get_prev_read_position), Some(previous));
    }

    // From the primary's first base to the read end: the supplementary's head is left.
    let mut model: AlignmentModel = original.clone();
    AlignmentModelEditor::new(&mut model).remove(459..1900);
    assert_eq!(model.get_read_sequence(), &sequence[..459]);
    assert_eq!(spans(&model), [(0, 458, Strand::Forward)]);
    assert_eq!(events(&model), [
        ((51, 52), AlignmentModelEventKind::Splicing),
        ((159, 160), AlignmentModelEventKind::Splicing),
        ((274, 275), AlignmentModelEventKind::Splicing)
    ]);
    assert!(model.get_bases().iter().all(|base| base.get_deletion_read_position().is_none()));
    for &(previous, next) in model.get_events().keys() {
        assert_eq!(model.get_event_at(previous).map(AlignmentModelEvent::get_next_read_position), Some(next));
        assert_eq!(model.get_event_at(next).map(AlignmentModelEvent::get_prev_read_position), Some(previous));
    }
}


/// Inserted bases take the read positions from the insertion point on. An event or record wholly
/// before that point stays, one wholly after it moves by the number of bases inserted, and one
/// that spans it widens. Inserted at the read start, the bases lie outside every record; appended
/// at the read end, they move nothing.
///
/// `scga-mini-rna-013` read 182: base 1000 lies inside both records, inside the breakpoint
/// (458, 1564), and between the deletion events (954, 955) and (1085, 1086).
#[test]
fn insert_unplaced_moves_what_follows_and_widens_what_spans_it() {
    let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_file.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-013-tumor_chunk_0000/182/ccs").unwrap();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
    let spans = |model: &AlignmentModel| -> Vec<(u32, u32, Strand)> {
        model.get_records().iter().map(|record| (record.read_start, record.read_end, record.reference_strand.clone())).collect()
    };
    let keys = |model: &AlignmentModel| -> Vec<(u32, u32)> {
        let mut keys: Vec<(u32, u32)> = model.get_events().keys().copied().collect();
        keys.sort();
        keys
    };
    assert_eq!(spans(&original), [(0, 1563, Strand::Forward), (459, 1893, Strand::Reverse)]);
    assert_eq!(keys(&original), [
        (51, 52), (159, 160), (274, 275), (458, 1564), (532, 533), (639, 640), (954, 955),
        (1085, 1086), (1252, 1253), (1472, 1473), (1515, 1516), (1656, 1657), (1850, 1851)
    ]);

    // Two bases in the middle of the read.
    let mut model: AlignmentModel = original.clone();
    AlignmentModelEditor::new(&mut model).insert_unplaced(1000, &[Nucleotide::G, Nucleotide::T], 3);
    assert_eq!(model.get_read_sequence(), format!("{}GT{}", &sequence[..1000], &sequence[1000..]));
    assert_eq!(spans(&model), [(0, 1565, Strand::Forward), (459, 1895, Strand::Reverse)]);
    assert_eq!(keys(&model), [
        (51, 52), (159, 160), (274, 275), (458, 1566), (532, 533), (639, 640), (954, 955),
        (1087, 1088), (1254, 1255), (1474, 1475), (1517, 1518), (1658, 1659), (1852, 1853)
    ]);
    assert_eq!(*model.get_event(458, 1566).unwrap().get_kind(), AlignmentModelEventKind::Breakpoint);
    for &(previous, next) in model.get_events().keys() {
        assert_eq!(model.get_event_at(previous).map(AlignmentModelEvent::get_next_read_position), Some(next));
        assert_eq!(model.get_event_at(next).map(AlignmentModelEvent::get_prev_read_position), Some(previous));
    }

    // One base at the read start: everything moves, and the new base belongs to no record.
    let mut model: AlignmentModel = original.clone();
    AlignmentModelEditor::new(&mut model).insert_unplaced(0, &[Nucleotide::G], 3);
    assert_eq!(model.get_read_sequence(), format!("G{sequence}"));
    assert_eq!(spans(&model), [(1, 1564, Strand::Forward), (460, 1894, Strand::Reverse)]);
    assert_eq!(keys(&model), keys(&original).into_iter().map(|(previous, next)| (previous + 1, next + 1)).collect::<Vec<(u32, u32)>>());

    // One base at the read end: nothing moves.
    let mut model: AlignmentModel = original.clone();
    AlignmentModelEditor::new(&mut model).insert_unplaced(1900, &[Nucleotide::G], 3);
    assert_eq!(model.get_read_sequence(), format!("{sequence}G"));
    assert_eq!(spans(&model), spans(&original));
    assert_eq!(keys(&model), keys(&original));
}


/// Restoring one deletion of a two-record read drops only that deletion's event. The restored
/// bases sit inside both records, which widen; every later event moves along, the flank moves
/// past the restored bases, and the other deletions stay parked on their flanks.
///
/// `scga-mini-rna-013` read 182: the Reverse primary parks a 3-base deletion at 1473 on flank
/// 1473 (event (1472, 1473)) and six 1-base deletions elsewhere.
#[test]
fn restore_deletion_drops_only_its_own_event() {
    let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_file.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-013-tumor_chunk_0000/182/ccs").unwrap();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
    let spans = |model: &AlignmentModel| -> Vec<(u32, u32, Strand)> {
        model.get_records().iter().map(|record| (record.read_start, record.read_end, record.reference_strand.clone())).collect()
    };
    let keys = |model: &AlignmentModel| -> Vec<(u32, u32)> {
        let mut keys: Vec<(u32, u32)> = model.get_events().keys().copied().collect();
        keys.sort();
        keys
    };
    let parked = |model: &AlignmentModel| -> Vec<u32> {
        model.get_bases().iter().filter_map(AlignmentModelBase::get_deletion_read_position).collect()
    };
    assert_eq!(parked(&original), [955, 1086, 1253, 1473, 1516, 1657, 1851]);
    assert_eq!(original.get_base(1473).get_deleted_reference_bases().len(), 3);
    assert_eq!(*original.get_event(1472, 1473).unwrap().get_kind(), AlignmentModelEventKind::Deletion);

    let mut model: AlignmentModel = original.clone();
    assert!(AlignmentModelEditor::new(&mut model).restore_deletion(1473, 7));
    assert_eq!(model.num_bases(), 1903);
    assert_eq!(spans(&model), [(0, 1566, Strand::Forward), (459, 1896, Strand::Reverse)]);
    assert_eq!(keys(&model), [
        (51, 52), (159, 160), (274, 275), (458, 1567), (532, 533), (639, 640), (954, 955),
        (1085, 1086), (1252, 1253), (1518, 1519), (1659, 1660), (1853, 1854)
    ]);
    assert_eq!(model.get_base(1476).get_placement(), original.get_base(1473).get_placement(), "the flank follows the restored bases");
    assert_eq!(parked(&model), [955, 1086, 1253, 1519, 1660, 1854]);
    for &(previous, next) in model.get_events().keys() {
        assert_eq!(model.get_event_at(previous).map(AlignmentModelEvent::get_next_read_position), Some(next));
        assert_eq!(model.get_event_at(next).map(AlignmentModelEvent::get_prev_read_position), Some(previous));
    }
}


/// A deletion parked between a kept flank and a cut has nothing left to restore it against, so
/// `retain` clears it. On the Reverse strand that flank is the first base kept and the deletion
/// is parked at the range start; on the Forward strand it is the last base kept and the deletion
/// is parked at the range end. Deletions parked inside the range stay.
///
/// `scga-mini-rna-013` read 182 parks seven deletions on its Reverse primary, the first at 955 on
/// flank 955. `scga-mini-rna-007` read 569 parks one at 766 on its Forward flank 765.
#[test]
fn retain_clears_a_deletion_parked_at_a_cut_on_the_flank_it_keeps() {
    for (bam_name, read_name, range, flank, parked_after) in [
        ("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam", "scga-mini-rna-013-tumor_chunk_0000/182/ccs", 955..1900u32, 955u32, vec![131u32, 298, 518, 561, 702, 896]),
        ("alignment/scga-mini-rna-007-tumor_minimap2_sorted.bam", "scga-mini-rna-007-tumor_chunk_0000/569/ccs", 0..766, 765, vec![])
    ] {
        let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(bam_name)).unwrap();
        let bam_file: &str = bam_file.to_str().unwrap();
        let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
        let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
        let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
        reader.read_header().unwrap();
        let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
        let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
        let deletion: u32 = original.get_base(flank).get_deletion_read_position().unwrap();
        assert!(deletion == range.start || deletion == range.end, "the deletion is parked at a cut");

        let mut model: AlignmentModel = original.clone();
        AlignmentModelEditor::new(&mut model).retain(range.clone());
        assert_eq!(model.get_read_sequence(), &sequence[range.start as usize..range.end as usize]);
        let kept: &AlignmentModelBase = model.get_base(flank - range.start);
        assert_eq!(kept.get_placement(), original.get_base(flank).get_placement(), "the flank is kept");
        assert!(kept.get_deleted_reference_bases().is_empty(), "its deletion is not");
        assert_eq!(model.get_bases().iter().filter_map(AlignmentModelBase::get_deletion_read_position).collect::<Vec<u32>>(), parked_after);
        for &(previous, next) in model.get_events().keys() {
            assert_eq!(model.get_event_at(previous).map(AlignmentModelEvent::get_next_read_position), Some(next));
            assert_eq!(model.get_event_at(next).map(AlignmentModelEvent::get_prev_read_position), Some(previous));
        }
    }
}


/// Edits with nothing to do leave every part of the model as it was: an empty removal, an empty
/// insertion, keeping the whole read, a restore where no deletion is parked, and a substitution
/// by the nucleotide the base already has, even at another quality.
///
/// `scga-mini-rna-013` read 182: base 1000 lies inside both records and parks no deletion.
#[test]
fn edits_with_nothing_to_do_leave_the_model_as_it_was() {
    let bam_file = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_file.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 1);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-013-tumor_chunk_0000/182/ccs").unwrap();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(&mut reader, read_id, &record_positions_map);
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let qualities: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let original: AlignmentModel = AlignmentModel::new(read_id, &sequence, &qualities, &records);
    let mut model: AlignmentModel = original.clone();
    let before: String = format!("{model:?}");
    let nucleotide: Nucleotide = model.get_base(1000).get_nucleotide().clone();
    assert!(model.get_bases().iter().all(|base| base.get_deletion_read_position() != Some(1000)));

    AlignmentModelEditor::new(&mut model).remove(1000..1000);
    assert_eq!(format!("{model:?}"), before, "an empty removal");
    AlignmentModelEditor::new(&mut model).insert_unplaced(1000, &[], 3);
    assert_eq!(format!("{model:?}"), before, "an empty insertion");
    AlignmentModelEditor::new(&mut model).retain(0..1900);
    assert_eq!(format!("{model:?}"), before, "the whole read kept");
    assert!(!AlignmentModelEditor::new(&mut model).restore_deletion(1000, 7));
    assert_eq!(format!("{model:?}"), before, "no deletion parked at 1000");
    assert!(!AlignmentModelEditor::new(&mut model).substitute(1000, nucleotide, 5));
    assert_eq!(format!("{model:?}"), before, "the same nucleotide at another quality");
}
