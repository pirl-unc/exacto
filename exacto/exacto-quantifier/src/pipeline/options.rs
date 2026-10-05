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


#[derive(Clone,Debug)]
pub struct QuantifyRNAAbundancesOptions {
    pub pseudo_count: f64,
    pub max_iter: usize,
    pub tol: f64
}

impl Default for QuantifyRNAAbundancesOptions {
    fn default() -> Self {
        Self {
            pseudo_count: 1e-2,
            max_iter: 10000,
            tol: 1e-9
        }
    }
}