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


use lazy_static::lazy_static;
use std::collections::HashMap;


pub const SINGLE_NUCLEOTIDE_VARIANT: &str = "SNV";
pub const MULTI_NUCLEOTIDE_VARIANT: &str = "MNV";
pub const INSERTION: &str = "INS";
pub const DELETION: &str = "DEL";
pub const INVERSION: &str = "INV";
pub const DUPLICATION: &str = "DUP";
pub const TRANSLOCATION: &str = "TRA";
pub const BREAKPOINT: &str = "BND";
pub const REFERENCE: &str = "REF";

pub const QUANTIFIER_ALL: &str = "all";
pub const QUANTIFIER_ANY: &str = "any";
pub const QUANTIFIER_MEDIAN: &str = "median";
pub const QUANTIFIER_AVERAGE: &str = "average";
pub const QUANTIFIER_MIN: &str = "min";
pub const QUANTIFIER_MAX: &str = "max";

pub const OPERATOR_LESS_THAN: &str = "<";
pub const OPERATOR_LESS_THAN_EQUAL_TO: &str = "<=";
pub const OPERATOR_GREATER_THAN: &str = ">";
pub const OPERATOR_GREATER_THAN_EQUAL_TO: &str = ">=";
pub const OPERATOR_EQUAL_TO: &str = "==";
pub const OPERATOR_NOT_EQUAL_TO: &str = "!=";
pub const OPERATOR_IN: &str = "in";

lazy_static! {
    pub static ref VARIANT_TYPES_MAP: HashMap<&'static str, String> = {
        let mut map = HashMap::new();
        map.insert(SINGLE_NUCLEOTIDE_VARIANT, SINGLE_NUCLEOTIDE_VARIANT.to_string());
        map.insert(MULTI_NUCLEOTIDE_VARIANT, MULTI_NUCLEOTIDE_VARIANT.to_string());
        map.insert(INSERTION, DUPLICATION.to_string() + ";" + INSERTION);
        map.insert(DUPLICATION, DUPLICATION.to_string() + ";" + INSERTION);
        map.insert(DELETION, DELETION.to_string());
        map.insert(INVERSION, BREAKPOINT.to_string() + ";" + INVERSION + ";" + TRANSLOCATION);
        map.insert(TRANSLOCATION, BREAKPOINT.to_string() + ";" + INVERSION + ";" + TRANSLOCATION);
        map.insert(BREAKPOINT, BREAKPOINT.to_string() + ";" + INVERSION + ";" + TRANSLOCATION);
        return map;
    };
}