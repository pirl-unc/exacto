from exactolib.constants import TranslationStrategy

from .data import get_data_path
from exactolib.main import translate_fastq_file


def test_translate_1():
    df_peptides = translate_fastq_file(
        fastq_file=get_data_path(name='fastq/rna-100-tumor_long-read.fastq.gz'),
        temp_dir='',
        strategy=TranslationStrategy.LONGEST_ORF,
        num_threads=1
    )
    assert(len(df_peptides) == 2)

