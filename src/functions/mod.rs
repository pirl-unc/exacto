pub mod annotate_variant_calls;
pub mod build_genome_variation_graph;
pub mod build_transcriptome_variation_graph;
pub mod cluster_rna_reads;
pub mod correct_rna_reads;
pub mod determine_rna_consensus;
pub mod identify_somatic_dna_variants;
pub mod identify_germline_dna_variants;
pub mod identify_rna_transcript_variants;
pub mod integrate_dna_rna_variants;
pub mod quantify_rna_abundances;
pub mod remove_unspliced_rnas;
pub mod stitch_reference_transcripts;
pub mod translate_fastx_file;
pub mod translate_sequence;
pub mod translate_transcripts;


/// The Python-facing error for a failed output write, so a wrapper can `?` the core writers'
/// results instead of discarding them.
pub(crate) fn io_error(error: impl std::fmt::Display) -> pyo3::PyErr {
    pyo3::PyErr::new::<pyo3::exceptions::PyIOError, _>(error.to_string())
}


/// The Python-facing error for a translation that stopped on its input: an unreadable or
/// unwritable file is an IOError, anything else about the input a ValueError.
pub(crate) fn translator_error(error: exacto::translator::prelude::TranslatorError) -> pyo3::PyErr {
    use exacto::translator::prelude::TranslatorError;
    match error {
        TranslatorError::File { .. } => pyo3::PyErr::new::<pyo3::exceptions::PyIOError, _>(error.to_string()),
        _ => pyo3::PyErr::new::<pyo3::exceptions::PyValueError, _>(error.to_string())
    }
}


/// The Python-facing error for a variation graph that could not be built: an unreadable file is an
/// IOError, a problem with the input table a ValueError.
pub(crate) fn graph_error(error: exacto::graph::prelude::GraphError) -> pyo3::PyErr {
    use exacto::graph::prelude::GraphError;
    match error {
        GraphError::File { .. } => pyo3::PyErr::new::<pyo3::exceptions::PyIOError, _>(error.to_string()),
        _ => pyo3::PyErr::new::<pyo3::exceptions::PyValueError, _>(error.to_string())
    }
}


/// The temp directory for the RNA pipelines, which hand `temp_dir` straight to
/// `NamedTempFile::new_in`, where an empty path is the current directory. An empty path means
/// TMPDIR, falling back to the system temp directory, as the DNA pipelines resolve it.
pub(crate) fn resolve_temp_dir(temp_dir: &str) -> String {
    if temp_dir.is_empty() {
        std::env::var("TMPDIR").unwrap_or_else(|_| std::env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    }
}
