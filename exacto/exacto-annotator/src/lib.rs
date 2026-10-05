extern crate bio;
extern crate exacto_caller;
extern crate exacto_core;
extern crate flate2;
extern crate indicatif;
extern crate log;
extern crate noodles_bgzf;
extern crate noodles_fasta;
extern crate once_cell;
extern crate phf;
extern crate polars;
extern crate rayon;
extern crate serde;
extern crate tempfile;

#[cfg(test)]
mod tests;
pub mod annotation;
pub mod pipeline;
pub mod prelude;