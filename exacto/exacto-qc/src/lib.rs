extern crate bimap;
extern crate exacto_caller;
extern crate exacto_core;
extern crate noodles_bam;
extern crate noodles_sam;
extern crate rayon;
extern crate serde;
extern crate tempfile;

#[cfg(test)]
mod tests;
pub mod pipeline;
pub mod io;
pub mod prelude;
pub mod common;