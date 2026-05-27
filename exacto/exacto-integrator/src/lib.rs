extern crate bimap;
extern crate csv;
extern crate exacto_annotator;
extern crate exacto_caller;
extern crate exacto_core;
extern crate indicatif;
extern crate polars;
extern crate rayon;

#[cfg(test)]
mod tests;
pub mod algorithms;
pub mod structs;
pub mod prelude;
pub mod common;
pub mod io;