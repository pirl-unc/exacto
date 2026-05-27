extern crate bimap;
extern crate bio;
extern crate bitvec;
extern crate chrono;
extern crate csv;
extern crate env_logger;
extern crate flate2;
extern crate log;
extern crate noodles_bam;
extern crate noodles_bgzf;
extern crate noodles_fasta;
extern crate noodles_util;
extern crate noodles_sam;
extern crate once_cell;
extern crate phf;
extern crate polars;
extern crate rayon;
extern crate serde;
extern crate sysinfo;
extern crate tempfile;

#[cfg(test)]
mod tests;
pub mod common;
pub mod index;
pub mod annotation;
pub mod prelude;
pub mod macros;
mod io;
