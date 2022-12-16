use pyo3::prelude::*;

mod variant_calling;


/// Identifies DNA variants in a BAM file.
#[pyfunction]
fn identify_dna_variants(bam_file: &str) -> PyResult<String> {
    Ok((1+2).to_string())
}


/// A Python module implemented in Rust.
#[pymodule]
fn exactors(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(identify_dna_variants, m)?)?;
    Ok(())
}