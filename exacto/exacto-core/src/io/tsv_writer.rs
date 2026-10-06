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


use csv::{Error, WriterBuilder};
use serde::Serialize;
use std::path::Path;


pub fn write_tsv_file<T, I>(records: I, path: &Path) -> Result<(), Error>
where
    T: Serialize,
    I: IntoIterator<Item = T>
{
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)?;
    for record in records {
        writer.serialize(record)?;
    }
    writer.flush()?;
    Ok(())
}


/// Write `records` as a tab-separated table, like `write_tsv_file`. A table without records is
/// written as its header line, the field names of `T`, so that any TSV reader reads it back as an
/// empty table: csv writes the header with the first record, so a default record gives it here.
pub fn write_tsv_table<T, I>(records: I, path: &Path) -> Result<(), Error>
where
    T: Serialize + Default,
    I: IntoIterator<Item = T>
{
    let mut records = records.into_iter().peekable();
    if records.peek().is_some() {
        return write_tsv_file(records, path);
    }
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(Vec::new());
    writer.serialize(T::default())?;
    writer.flush()?;
    let table: &[u8] = writer.get_ref();
    let header_end: usize = table.iter().position(|&byte| byte == b'\n').map_or(table.len(), |index| index + 1);
    std::fs::write(path, &table[..header_end])?;
    Ok(())
}
