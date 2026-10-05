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


use thiserror::Error;


/// Why a graph could not be built. Each variant names the file, column, row or transcript at fault.
#[derive(Debug, Error)]
pub enum GraphError {
    #[error("{file} could not be read: {reason}")]
    File { file: Box<str>, reason: Box<str> },

    #[error("column {column}: {reason}")]
    Column { column: Box<str>, reason: Box<str> },

    #[error("{row}: {reason}")]
    Row { row: Box<str>, reason: Box<str> },

    #[error("transcript {name}: {reason}")]
    Transcript { name: Box<str>, reason: Box<str> }
}
