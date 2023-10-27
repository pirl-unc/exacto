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
use pyo3::prelude::*;
use serde::{Serialize, Deserialize};
use serde_json::Value;
use crate::constants;
use crate::variant_call::VariantCall;
use crate::variant::Variant;


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantFilter {
    pub quantifier: String,
    pub attribute: String,
    pub operator: String,
    pub value: Value,
    pub sample_ids: Vec<String>
}

impl VariantFilter {
    pub fn new(
        quantifier: String,
        attribute: String,
        operator: String,
        value: Value,
        sample_ids: Vec<String>) -> Self {
        Self {
            quantifier: quantifier,
            attribute: attribute,
            operator: operator,
            value: value,
            sample_ids: sample_ids
        }
    }

    pub fn keep_variant_call(&self, variant_call: &VariantCall) -> bool {
        // Check if VariantCall is eligible for filtering by this VariantFilter.
        // Return true if VariantCall is ineligible for filtering.
        if self.sample_ids.contains(&variant_call.sample_id) == false {
            return true;
        }

        let string_attributes = vec![
            "id",
            "source_id",
            "sample_id",
            "phase_block_id",
            "clone_id",
            "nucleic_acid",
            "variant_calling_method",
            "sequencing_platform",
            "chromosome_1",
            "chromosome_2",
            "reference_allele",
            "alternate_allele",
            "filter",
            "variant_type",
            "variant_subtype"
        ];
        let numeric_attributes = vec![
            "position_1",
            "position_2",
            "quality_score",
            "variant_size",
            "total_read_count",
            "reference_allele_read_count",
            "alternate_allele_read_count",
            "alternate_allele_fraction"
        ];
        let boolean_attributes = vec![
            "precise"
        ];
        if string_attributes.contains(&self.attribute.as_str()) {
            let attribute_value: &str = match self.attribute.as_str() {
                "id" => &variant_call.id.as_str(),
                "source_id" => variant_call.source_id.as_str(),
                "sample_id" => variant_call.sample_id.as_str(),
                "phase_block_id" => variant_call.phase_block_id.as_str(),
                "clone_id" => variant_call.clone_id.as_str(),
                "nucleic_acid" => variant_call.nucleic_acid.as_str(),
                "variant_calling_method" => variant_call.variant_calling_method.as_str(),
                "sequencing_platform" => variant_call.sequencing_platform.as_str(),
                "chromosome_1" => variant_call.chromosome_1.as_str(),
                "chromosome_2" => variant_call.chromosome_2.as_str(),
                "reference_allele" => variant_call.reference_allele.as_str(),
                "alternate_allele" => variant_call.alternate_allele.as_str(),
                "filter" => variant_call.filter.as_str(),
                "variant_type" => variant_call.variant_type.as_str(),
                "variant_subtype" => variant_call.variant_subtype.as_str(),
                _ => "",
            };
            match &self.value {
                Value::String(filter_value) => {
                    if self.operator == constants::OPERATOR_EQUAL_TO {
                        if attribute_value == filter_value {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_NOT_EQUAL_TO {
                        if attribute_value != filter_value {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Array(filter_value) => {
                    if self.operator == constants::OPERATOR_IN {
                        let attribute_value_as_value = Value::String(attribute_value.to_string());
                        if filter_value.contains(&attribute_value_as_value) {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string array: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Bool(filter_value) => {
                    eprintln!("Filter value cannot be a boolean for a string attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Number(filter_value) => {
                    eprintln!("Filter value cannot be a number for a string attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Null => {
                    eprintln!("Filter value cannot be null for a string attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Object(_) => {
                    eprintln!("Filter value cannot be an object for a string attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
            }
        } else if numeric_attributes.contains(&self.attribute.as_str()) {
            let attribute_value: f64 = match self.attribute.as_str() {
                "position_1" => variant_call.position_1 as f64,
                "position_2" => variant_call.position_2 as f64,
                "quality_score" => variant_call.quality_score as f64,
                "variant_size" => variant_call.variant_size as f64,
                "total_read_count" => variant_call.total_read_count as f64,
                "reference_allele_read_count" => variant_call.reference_allele_read_count as f64,
                "alternate_allele_read_count" => variant_call.alternate_allele_read_count as f64,
                "alternate_allele_fraction" => variant_call.alternate_allele_fraction as f64,
                _ => -1.0,
            };
            match &self.value {
                Value::String(filter_value) => {
                    let filter_value_: Result<f64, std::num::ParseFloatError> = filter_value.as_str().parse();
                    let filter_value_f64: f64 = match filter_value_ {
                        Ok(value) => value,
                        Err(_) => {
                            eprintln!("{}", format!("Cannot parse value to a f64 value: {}", self.value));
                            std::process::exit(exitcode::DATAERR);
                        },
                    };
                    if self.operator == constants::OPERATOR_EQUAL_TO {
                        if attribute_value == filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_LESS_THAN {
                        if attribute_value < filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_LESS_THAN_EQUAL_TO {
                        if attribute_value <= filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_GREATER_THAN {
                        if attribute_value > filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_GREATER_THAN_EQUAL_TO {
                        if attribute_value >= filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_NOT_EQUAL_TO {
                        if attribute_value != filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Array(filter_value) => {
                    eprintln!("Filter value cannot be a string array for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Bool(filter_value) => {
                    eprintln!("Filter value cannot be a boolean for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Number(filter_value) => {
                    let filter_value_f64 = filter_value.as_f64().unwrap();
                    if self.operator == constants::OPERATOR_EQUAL_TO {
                        if attribute_value == filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_LESS_THAN {
                        if attribute_value < filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_LESS_THAN_EQUAL_TO {
                        if attribute_value <= filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_GREATER_THAN {
                        if attribute_value > filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_GREATER_THAN_EQUAL_TO {
                        if attribute_value >= filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_NOT_EQUAL_TO {
                        if attribute_value != filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Null => {
                    eprintln!("Filter value cannot be null for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Object(_) => {
                    eprintln!("Filter value cannot be an object for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
            }
        } else if boolean_attributes.contains(&self.attribute.as_str()) {
             let attribute_value: bool = match self.attribute.as_str() {
                "precise" => variant_call.precise,
                _ => false,
            };
            match &self.value {
                Value::String(filter_value) => {
                    let filter_value_: Result<bool, std::str::ParseBoolError> = filter_value.as_str().parse();
                    let filter_value_bool: bool = match filter_value_ {
                        Ok(value) => value,
                        Err(_) => {
                            eprintln!("{}", format!("Cannot parse value to a boolean value: {}", self.value));
                            std::process::exit(exitcode::DATAERR);
                        },
                    };
                    if self.operator == constants::OPERATOR_EQUAL_TO {
                        if attribute_value == filter_value_bool {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_NOT_EQUAL_TO {
                        if attribute_value != filter_value_bool {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Array(filter_value) => {
                    eprintln!("Filter value cannot be a string array for a boolean attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Bool(filter_value) => {
                    let filter_value: bool = filter_value.clone();
                    if self.operator == constants::OPERATOR_EQUAL_TO {
                        if attribute_value == filter_value {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_NOT_EQUAL_TO {
                        if attribute_value != filter_value {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Number(filter_value) => {
                    eprintln!("Filter value cannot be a number for a boolean attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Null => {
                    eprintln!("Filter value cannot be null for a boolean attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Object(_) => {
                    eprintln!("Filter value cannot be an object for a boolean attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
            }
        }
        else {
            return false;
        }
    }

    pub fn keep(&self, variant: &Variant) -> bool {
        if self.quantifier == constants::QUANTIFIER_ALL {
            for variant_call in &variant.variant_calls {
                if self.keep_variant_call(&variant_call) == false {
                    return false;
                }
            }
            return true;
        } else if self.quantifier == constants::QUANTIFIER_ANY {
            for variant_call in &variant.variant_calls {
                if self.keep_variant_call(variant_call) {
                    return true;
                }
            }
            return false;
        } else if self.quantifier == constants::QUANTIFIER_AVERAGE ||
            self.quantifier == constants::QUANTIFIER_MEDIAN ||
            self.quantifier == constants::QUANTIFIER_MIN ||
            self.quantifier == constants::QUANTIFIER_MAX {
            let attribute_value: f64 = match self.quantifier.as_str() {
                constants::QUANTIFIER_AVERAGE => variant.get_mean(self.attribute.to_string()),
                constants::QUANTIFIER_MEDIAN => variant.get_median(self.attribute.to_string()),
                constants::QUANTIFIER_MIN => variant.get_min(self.attribute.to_string()),
                constants::QUANTIFIER_MAX => variant.get_max(self.attribute.to_string()),
                _ => -1.0,
            };
            match &self.value {
                Value::String(filter_value) => {
                    eprintln!("Filter value cannot be a string for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Array(filter_value) => {
                    eprintln!("Filter value cannot be a string array for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Bool(filter_value) => {
                    eprintln!("Filter value cannot be a boolean for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Number(filter_value) => {
                    let filter_value_f64 = filter_value.as_f64().unwrap();
                    if self.operator == constants::OPERATOR_EQUAL_TO {
                        if attribute_value == filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_LESS_THAN {
                        if attribute_value < filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_LESS_THAN_EQUAL_TO {
                        if attribute_value <= filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_GREATER_THAN {
                        if attribute_value > filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_GREATER_THAN_EQUAL_TO {
                        if attribute_value >= filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else if self.operator == constants::OPERATOR_NOT_EQUAL_TO {
                        if attribute_value != filter_value_f64 {
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        eprintln!("{}", format!("Unsupported operator for string: {}", self.operator));
                        std::process::exit(exitcode::DATAERR);
                    }
                }
                Value::Null => {
                    eprintln!("Filter value cannot be null for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
                Value::Object(_) => {
                    eprintln!("Filter value cannot be an object for a numeric attribute.");
                    std::process::exit(exitcode::DATAERR);
                }
            }
        } else {
            eprintln!("{}", format!("Unsupported quantifier: {}", self.quantifier));
            std::process::exit(exitcode::DATAERR);
        }
    }
}

impl Clone for VariantFilter {
    fn clone(&self) -> Self {
        let mut sample_ids: Vec<String> = Vec::new();
        for sample_id in &self.sample_ids {
            sample_ids.push(sample_id.clone());
        }
        VariantFilter {
            quantifier: self.quantifier.to_string(),
            attribute: self.attribute.to_string(),
            operator: self.operator.to_string(),
            value: self.value.clone(),
            sample_ids: sample_ids
        }
    }
}
