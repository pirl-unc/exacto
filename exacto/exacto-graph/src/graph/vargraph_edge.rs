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


use std::any::Any;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraphEdge {
    pub orientation_from: VarGraphOrientations,
    pub orientation_to: VarGraphOrientations,
    pub enabled: bool,
    pub annotations: HashMap<Box<str>, Arc<dyn Any + Send + Sync + 'static>>
}

impl PartialEq for VarGraphEdge {
    fn eq(&self, other: &Self) -> bool {
        self.orientation_from == other.orientation_from && 
            self.orientation_to == other.orientation_to &&
            self.enabled == other.enabled
    }
}

impl Hash for VarGraphEdge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.orientation_from.hash(state);
        self.orientation_to.hash(state);
        self.enabled.hash(state);
    }
}

impl VarGraphEdge {
    pub fn new(
        orientation_from: VarGraphOrientations,
        orientation_to: VarGraphOrientations
    ) -> Self {
        VarGraphEdge {
            orientation_from: orientation_from,
            orientation_to: orientation_to,
            enabled: true,
            annotations: HashMap::new()
        }
    }

    pub fn add_annotation(&mut self, key: &str, value: Arc<dyn Any + Send + Sync + 'static>) {
        self.annotations.insert(key.into(), value);
    }

    pub fn enable(&mut self) {
        self.enabled = true;
    }

    pub fn disable(&mut self) {
        self.enabled = false;
    }

    pub fn is_enabled(&self) -> bool {
        self.enabled
    }
}

impl Clone for VarGraphEdge {
    fn clone(&self) -> VarGraphEdge {
        VarGraphEdge {
            orientation_from: self.orientation_from.clone(),
            orientation_to: self.orientation_to.clone(),
            enabled: self.enabled,
            annotations: self.annotations.clone(),
        }
    }
}
