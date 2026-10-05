use super::*;


#[test]
fn to_string_formats_single_gene_transcript_exon_1() {
    let mut position_annotation: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    position_annotation.add_reference_gene_id("ENSG001".into());
    position_annotation.add_reference_transcript_id("ENSG001".into(), "ENST001".into());
    position_annotation.add_reference_exon_id("ENST001".into(), "EXON001".into());

    assert_eq!(position_annotation.to_string(), "ENSG001|ENSG001-ENST001|ENST001-EXON001");
}

#[test]
fn to_string_formats_single_gene_transcript_exon_2() {
    let mut position_annotation: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    position_annotation.add_reference_gene_id("ENSG001".into());
    position_annotation.add_reference_transcript_id("ENSG001".into(), "ENST001".into());
    position_annotation.add_reference_transcript_id("ENSG001".into(), "ENST002".into());
    position_annotation.add_reference_exon_id("ENST001".into(), "EXON001".into());
    position_annotation.add_reference_exon_id("ENST002".into(), "EXON002".into());

    assert_eq!(position_annotation.to_string(), "ENSG001|ENSG001-ENST001;ENSG001-ENST002|ENST001-EXON001;ENST002-EXON002");
}

#[test]
fn to_string_formats_single_gene_transcript_exon_3() {
    let mut position_annotation_1: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    position_annotation_1.add_reference_gene_id("ENSG001".into());
    position_annotation_1.add_reference_transcript_id("ENSG001".into(), "ENST001".into());
    position_annotation_1.add_reference_transcript_id("ENSG001".into(), "ENST002".into());
    position_annotation_1.add_reference_exon_id("ENST001".into(), "EXON001".into());
    position_annotation_1.add_reference_exon_id("ENST002".into(), "EXON002".into());

    let position_annotation_2: PositionAnnotation = position_annotation_1.clone();

    assert_eq!(position_annotation_2, position_annotation_1);
}

#[test]
fn to_string_formats_single_gene_transcript_exon_4() {
    let position_annotation_1: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    assert_eq!(position_annotation_1.to_string(), "||");

    let position_annotation_2: PositionAnnotation = PositionAnnotation::from_string("ENSG001|ENSG001-ENST001|");
    assert_eq!(position_annotation_2.genic_region, GenicRegion::Intronic);

    let position_annotation_3: PositionAnnotation = PositionAnnotation::from_string("||");
    assert_eq!(position_annotation_3.genic_region, GenicRegion::Intergenic);
}