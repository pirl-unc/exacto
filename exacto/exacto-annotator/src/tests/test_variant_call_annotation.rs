use exacto_caller::prelude::VariantType;
use exacto_core::prelude::*;
use polars::prelude::*;

use crate::prelude::*;


#[test]
fn test_variant_call_annotation_1() {
    let variant_call_annotation_1: VariantCallAnnotation = VariantCallAnnotation::new(
        1,
        "chr1",
        1001,
        "chr1",
        1003,
        VariantType::SingleNucleotideVariant,
        "A",
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic)
    );

    let variant_call_annotation_2: VariantCallAnnotation = variant_call_annotation_1.clone();

    assert_eq!(variant_call_annotation_1, variant_call_annotation_2);
}
