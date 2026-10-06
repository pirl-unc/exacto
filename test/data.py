import gzip
import os
from typing import Dict


def get_data_path(name: str) -> str:
    """
    Return the absolute path to a file in the test/data directory.

    Parameters
    ----------
        name    :   Name of file.

    Returns
    -------
        Absolute path to a file in the test/data directory.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)


def get_simulated_read_origins(sample_id: str) -> Dict[str, str]:
    """
    Return the simulated transcript each read of an scga-mini RNA sample was drawn from.

    Parameters
    ----------
        sample_id   :   Sample ID (e.g. 'scga-mini-rna-001-tumor').

    Returns
    -------
        Dictionary of read name (as in the aligned BAM, '<chunk>/<zmw>/ccs') to the name of
        the row of simulation/transcript/<sample_id>.transcript it was simulated from.
    """
    read_origins = {}
    for chunk in ('0000', '0001'):
        maf_file = get_data_path(name='simulation/fastq/pbsim3/%s/pbsim3/%s_chunk_%s_pbsim3.maf.gz' % (sample_id, sample_id, chunk))
        with gzip.open(maf_file, 'rt') as handle:
            # Each alignment block holds the transcript's 's' line, then the read's ('<chunk>/<zmw>/<n>').
            names = [line.split()[1] for line in handle if line.startswith('s ')]
        for transcript_name, read_name in zip(names[0::2], names[1::2]):
            read_origins[read_name.rsplit('/', 1)[0] + '/ccs'] = transcript_name
    return read_origins
