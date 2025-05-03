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


use serde::{Serialize, Deserialize};

use crate::prelude::IntegrationRecord;


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct IntegrationRecordSet {
    pub integration_records: HashSet<IntegrationRecord>

}

impl Hash for IntegrationRecordSet {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.integration_records.hash(state);
    }
}

impl IntegrationRecordSet {
    pub fn new() -> Self {
        Self {
            integration_records: HashSet::new()
        }
    }

    pub fn add_integration_record(&mut self, integration_record: IntegrationRecord) {
        self.integration_records.insert(integration_record);
    }
}

impl Clone for IntegrationRecordSet {
    fn clone(&self) -> Self {
        IntegrationRecord {
            integration_records: self.integration_records.clone()
        }
    }
}
