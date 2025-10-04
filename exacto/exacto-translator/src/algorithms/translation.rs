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


use rayon::prelude::*;

use crate::prelude::*;


pub fn translate_rnas(
    rnas: Vec<RNA>,
    num_threads: usize
) -> TranslationSet {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let translations: Vec<Translation> = thread_pool.install(|| {
        rnas
            .par_iter()
            .filter_map(|rna| rna.translate()) // Filters out `None` and unwraps `Some`.
            .collect()
    });
    let mut translation_set: TranslationSet = TranslationSet::new();
    for translation in translations {
        translation_set.add_translation(translation);
    }
    translation_set
}
