use exacto_core::prelude::GenicRegion;

use super::*;


#[test]
fn variant_call_annotation_clone_equals_original() {
    let variant_call_annotation_1: VariantCallAnnotation = VariantCallAnnotation::new(
        1,
        "chr1".into(),
        1001,
        "chr1".into(),
        1003,
        VariantType::SingleNucleotideVariant,
        "A".into(),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic)
    );

    let variant_call_annotation_2: VariantCallAnnotation = variant_call_annotation_1.clone();

    assert_eq!(variant_call_annotation_1, variant_call_annotation_2);
}
