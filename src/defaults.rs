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


// Minimum insertion size proportion between two insertions
// size proportion = smaller insertion size / longer insertion size
pub const MIN_INS_SIZE_PROPORTION: f64 = 0.5;

// Maximum insertion normalized edit (levenshtein) distance
// normalized edit distance = edit distance / longer insertion size
pub const MAX_INS_NORM_EDIT_DISTANCE: f64 = 0.5;

// Minimum deletion size proportion between two deletions
// size proportion = smaller deletion / longer deletion
pub const MIN_DEL_SIZE_PROPORTION: f64 = 0.5;
