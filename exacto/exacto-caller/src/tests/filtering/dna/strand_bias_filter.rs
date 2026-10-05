use super::*;


#[test]
fn has_strand_bias_returns_matches() {
    let strand_bias_exists: bool = has_strand_bias(
        0,
        8,
        12,
        13,
        0.05
    );
    assert_eq!(strand_bias_exists, true);
}
