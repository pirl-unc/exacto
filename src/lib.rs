// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


extern crate exitcode;
extern crate pyo3;
extern crate serde_json;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
mod constants;
mod genomic_range;
mod genomic_ranges_list;
mod utilities;
mod variant;
mod variant_annotation;
mod variant_call;
mod variant_filter;
mod variants_list;
use genomic_range::GenomicRange;
use genomic_ranges_list::GenomicRangesList;
use variant::Variant;
use variant_annotation::VariantAnnotation;
use variant_call::VariantCall;
use variant_filter::VariantFilter;
use variants_list::VariantsList;


/// This function deserializes a list of serialized VariantsList objects.
///
/// # Arguments
/// * `py_list`         -   list of serialized VariantsList objects.
///
/// # Returns
/// * `variants_list`   -   vector of VariantsList objects.
fn deserialize_variants_lists(py_list: &PyList) -> Vec<VariantsList> {
    let mut variants_lists: Vec<VariantsList> = Vec::new();
    for py_str in py_list.iter() {
        variants_lists.push(deserialize_variants_list(&py_str.to_string()));
    }
    return variants_lists;
}

/// This function deserializes a VariantsList object.
///
/// # Arguments
/// * `json_str`        -   serialized VariantsList object.
///
/// # Returns
/// * `variants_list`   -   VariantsList object.
fn deserialize_variants_list(json_str: &str) -> VariantsList {
    let variants_list: Result<VariantsList, serde_json::Error> = serde_json::from_str(json_str);
    match variants_list {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            eprintln!("Error deserializing JSON: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    }
}

/// This function deserializes a list of VariantFilter objects.
///
/// # Arguments
/// * `py_list`             -   a list of serialized VariantsFilter objects.
///
/// # Returns
/// * `variant_filters`     -   a vector of VariantFilter objects.
fn deserialize_variant_filters(py_list: &PyList) -> Vec<VariantFilter> {
    let mut variant_filters: Vec<VariantFilter> = Vec::new();
    for py_str in py_list.iter() {
        variant_filters.push(deserialize_variant_filter(&py_str.to_string()));
    }
    return variant_filters;
}

fn deserialize_variant_filter(json_str: &str) -> VariantFilter {
    let variant_filter: Result<VariantFilter, serde_json::Error> = serde_json::from_str(json_str);
    match variant_filter {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            eprintln!("Error deserializing JSON: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    }
}

/// This function merges a vector of serialized VariantsList objects into one.
///
/// # Arguments
/// * `py_list`                 -   list of serialized VariantsList objects.
/// * `num_threads`             -   number of threads.
/// * `max_neighbor_distance`   -   maximum neighbor distance.
///
/// # Returns
/// * A serialized VariantsList object.
#[pyfunction]
fn merge_variants_lists(
    py: Python, py_list: &PyList,
    num_threads: usize,
    max_neighbor_distance: isize) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let mut variants_lists: Vec<VariantsList> = deserialize_variants_lists(py_list);

    // Step 2. Merge VariantsList objects
    let merged_variants_list: VariantsList = VariantsList::merge(
        variants_lists,
        num_threads,
        max_neighbor_distance,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize merged VariantsList object
    let serialized = serde_json::to_string(&merged_variants_list).expect("Serialization of the merged VariantsList object failed");

    Ok(serialized)
}

/// This function filters a serialized VariantsList object and returns a filtered VariantsList.
///
/// # Arguments
/// * `py_str`                  -   serialized VariantsList object.
/// * `py_list`                 -   a list of serialized VariantFilter objects.
/// * `num_threads`             -   number of threads.
///
/// # Returns
/// * A serialized VariantsList string.
#[pyfunction]
fn filter_variants_list(
    py: Python, py_str:
    String, py_list: &PyList,
    num_threads: usize) -> PyResult<String> {
    // Step 1. Deserialize VariantsList object
    let mut variants_list: VariantsList = deserialize_variants_list(&py_str);

    // Step 2. Deserialize VariantFilter objects
    let variant_filters: Vec<VariantFilter> = deserialize_variant_filters(py_list);

    // Step 3. Filter VariantsList object
    let filtered_variants_list: VariantsList = variants_list.filter(variant_filters, num_threads);

    // Step 4. Serialize filtered VariantsList object
    let serialized = serde_json::to_string(&filtered_variants_list).expect("Serialization of the filtered VariantsList object failed");

    Ok(serialized)
}

/// This function identifies nearby Variant objects.
///
/// # Arguments
/// * `target_variants_list_str`    -   serialized VariantsList object.
/// * `query_variants_list_str`     -   serialized VariantsList object.
/// * `num_threads`                 -   number of threads.
/// * `max_neighbor_distance`       -   maximum neighbor distance.
///
/// # Returns
/// * `nearby_variants_map`         -   HashMap where key is Variant.id and
///                                     value is a vector of query Variant.id
#[pyfunction]
fn find_nearby_variants(
    py: Python,
    target_variants_list_str: String,
    query_variants_list_str: String,
    num_threads: usize,
    max_neighbor_distance: isize) -> Py<PyAny> {
    // Step 1. Deserialize VariantsList objects
    let mut target_variants_list: VariantsList = deserialize_variants_list(&target_variants_list_str);
    let mut query_variants_list: VariantsList = deserialize_variants_list(&query_variants_list_str);

    // Step 2. Find nearby variants
    let nearby_variants_map: HashMap<String, Vec<String>> = target_variants_list.find_nearby_variants(
        query_variants_list,
        max_neighbor_distance,
        num_threads,
        &constants::VARIANT_TYPES_MAP
    );

    return Python::with_gil(|py| {
        nearby_variants_map.to_object(py)
    });
}

#[pymodule]
fn exactors(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(merge_variants_lists, m)?);
    m.add_function(wrap_pyfunction!(filter_variants_list, m)?);
    m.add_function(wrap_pyfunction!(find_nearby_variants, m)?);
    Ok(())
}
