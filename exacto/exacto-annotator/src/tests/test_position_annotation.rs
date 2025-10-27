use exacto_core::prelude::*;
use polars::prelude::*;

use crate::prelude::*;


#[test]
fn test_position_annotation_1() {
    let mut position_annotation: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    position_annotation.add_gene_id("ENSG001");
    position_annotation.add_transcript_id("ENSG001", "ENST001");
    position_annotation.add_exon_id("ENST001", "EXON001");

    assert_eq!(position_annotation.to_string(), "ENSG001;ENSG001-ENST001;ENST001-EXON001");
}

#[test]
fn test_position_annotation_2() {
    let mut position_annotation: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    position_annotation.add_gene_id("ENSG001");
    position_annotation.add_transcript_id("ENSG001", "ENST001");
    position_annotation.add_transcript_id("ENSG001", "ENST002");
    position_annotation.add_exon_id("ENST001", "EXON001");
    position_annotation.add_exon_id("ENST002", "EXON002");

    assert_eq!(position_annotation.to_string(), "ENSG001;ENSG001-ENST001,ENSG001-ENST002;ENST001-EXON001,ENST002-EXON002");
}

#[test]
fn test_position_annotation_3() {
    let mut position_annotation_1: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    position_annotation_1.add_gene_id("ENSG001");
    position_annotation_1.add_transcript_id("ENSG001", "ENST001");
    position_annotation_1.add_transcript_id("ENSG001", "ENST002");
    position_annotation_1.add_exon_id("ENST001", "EXON001");
    position_annotation_1.add_exon_id("ENST002", "EXON002");

    let position_annotation_2: PositionAnnotation = position_annotation_1.clone();

    assert_eq!(position_annotation_2, position_annotation_1);
}

#[test]
fn test_position_annotation_4() {
    let position_annotation_1: PositionAnnotation = PositionAnnotation::new(GenicRegion::Exonic);
    assert_eq!(position_annotation_1.to_string(), ";;");

    let position_annotation_2: PositionAnnotation = PositionAnnotation::from_string("ENSG001;ENSG001-ENST001;");
    assert_eq!(position_annotation_2.genic_region, GenicRegion::Intronic);

    let position_annotation_3: PositionAnnotation = PositionAnnotation::from_string(";;");
    assert_eq!(position_annotation_3.genic_region, GenicRegion::Intergenic);
}