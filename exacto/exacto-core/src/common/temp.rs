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


use std::env;
use std::path::PathBuf;


pub fn fetch_temp_directory(temp_dir: &str) -> PathBuf {
    let dir: PathBuf = if temp_dir.is_empty() {
        env::temp_dir()
    } else {
        PathBuf::from(temp_dir)
    };
    if !dir.is_dir() {
        panic!("Directory does not exist: {}", dir.display());
    }
    dir
}