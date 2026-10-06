import os

import pandas as pd
import pysam
import random
from typing import Dict, List, Tuple
from vstolib.gencode import Gencode


GENOME_FASTA_FILE = "/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz"
GENCODE_GTF_FILE = "/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/gencode.v41.annotation.chr17-18.gtf.gz"
FASTA_DIR = '../../../../test/data/simulation/fasta'
TRANSCRIPT_DIR = '../../../../test/data/simulation/transcript'
NUM_SENSE_READS = 200
NUM_ANTISENSE_READS = 0
NUM_PARTIALLY_SPLICED_TRANSCRIPTS = 3
NUM_SENSE_READS_PARTIALLY_SPLICED = 10
MIN_RETAINED_INTRONS = 1
MAX_RETAINED_INTRONS = 3
RANDOM_SEED = 1
REFERENCE_TRANSCRIPT_IDS: Dict[str, List[str]] = {
    'scga-mini-rna-001-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-002-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-003-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-004-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-005-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-006-tumor': ['ENST00000698742'],                                     # WRAP53
    'scga-mini-rna-007-tumor': ['ENST00000263080', 'ENST00000317744'],                  # ASPA, WSCD1
    'scga-mini-rna-008-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-009-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-010-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-011-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-012-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-013-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-014-tumor': ['ENST00000269305'],                                     # TP53
    'scga-mini-rna-015-tumor': ['ENST00000263080', 'ENST00000317744',
                                'ENST00000575425'],                                     # ASPA, WSCD1, ACAP1
    'scga-mini-rna-016-tumor': ['ENST00000263080', 'ENST00000317744'],                  # ASPA, WSCD1
}


def reverse_complement(sequence: str) -> str:
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'N': 'N',
        'n': 'n'
    }
    reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
    return reverse_complement_seq


def build_reference_transcript_dict(
        gencode: Gencode,
        transcript_ids_stable: List[str]
) -> Dict[str, Dict]:
    """
    Build a lookup of the reference transcripts to simulate.

    Exons are sorted by genomic start ascending rather than by exon number, so that intron i is
    always (exons[i].end + 1, exons[i+1].start - 1) no matter which strand the gene is on. The
    strand is applied once, at the end, when the sequence is assembled.

    Parameters:
        gencode                 :   Gencode object.
        transcript_ids_stable   :   List of stable Ensembl transcript IDs (no version suffix).

    Returns:
        Dictionary of stable transcript ID to {'chromosome', 'strand', 'exons'}.
    """
    reference_transcripts = {}
    for transcript_id_stable in transcript_ids_stable:
        df_transcript = gencode.df_transcripts[
            gencode.df_transcripts['transcript_id_stable'] == transcript_id_stable
        ]
        assert len(df_transcript) == 1, \
            '%s matched %i transcripts in %s' % (transcript_id_stable, len(df_transcript), GENCODE_GTF_FILE)
        df_exons = gencode.df_exons[
            gencode.df_exons['transcript_id'] == df_transcript['transcript_id'].values[0]
        ].sort_values(by=['start'], ascending=True)
        assert len(df_exons) >= 2, \
            '%s has %i exon(s) and so has no introns to retain' % (transcript_id_stable, len(df_exons))
        reference_transcripts[transcript_id_stable] = {
            'chromosome': str(df_exons['chromosome'].values[0]),
            'strand': str(df_exons['strand'].values[0]),
            'exons': [(int(start), int(end)) for start, end in zip(df_exons['start'], df_exons['end'])]
        }
    return reference_transcripts


def get_introns(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Return the introns between genomically sorted exons, in the same order.
    """
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]


def get_transcript_sequence(
        fasta: pysam.FastaFile,
        reference_transcript: Dict,
        retained_intron_indices: Tuple[int, ...]
) -> str:
    """
    Assemble a transcript sequence, retaining the introns at the given indices.

    An empty retained_intron_indices yields the fully spliced (mature) transcript.

    Parameters:
        fasta                   :   Reference genome FASTA file.
        reference_transcript    :   Entry from build_reference_transcript_dict.
        retained_intron_indices :   Indices into get_introns(exons) of the introns left unspliced.

    Returns:
        Transcript sequence, on the transcript's own strand.
    """
    chromosome = reference_transcript['chromosome']
    exons = reference_transcript['exons']
    introns = get_introns(exons)

    # GENCODE coordinates are 1-based inclusive; pysam fetch is 0-based half-open.
    blocks = []
    for i, (start, end) in enumerate(exons):
        blocks.append(str(fasta.fetch(chromosome, start - 1, end)))
        if i < len(introns) and i in retained_intron_indices:
            intron_start, intron_end = introns[i]
            blocks.append(str(fasta.fetch(chromosome, intron_start - 1, intron_end)))
    sequence = ''.join(blocks)

    if reference_transcript['strand'] == '-':
        sequence = reverse_complement(sequence)

    return sequence


def sample_retained_intron_indices(
        num_introns: int,
        num_samples: int,
        rng: random.Random
) -> List[Tuple[int, ...]]:
    """
    Draw distinct sets of introns to retain.

    Each set holds between MIN_ and MAX_RETAINED_INTRONS introns (capped at what the transcript
    has). Duplicate sets are rejected rather than accepted, so the isoforms of one transcript are
    always different from one another. Fewer than num_samples sets are returned only if the
    transcript does not have enough distinct combinations, which is checked by the caller.
    """
    max_retained = min(MAX_RETAINED_INTRONS, num_introns)
    combinations = set()
    ordered = []
    # Bounded so that a transcript with too few combinations terminates instead of spinning.
    for _ in range(num_samples * 1000):
        if len(ordered) == num_samples:
            break
        num_retained = rng.randint(MIN_RETAINED_INTRONS, max_retained)
        indices = tuple(sorted(rng.sample(range(num_introns), num_retained)))
        if indices in combinations:
            continue
        combinations.add(indices)
        ordered.append(indices)
    return ordered


def create_partially_spliced_transcripts(
        fasta: pysam.FastaFile,
        reference_transcripts: Dict[str, Dict],
        sample_id: str,
        transcript_ids_stable: List[str]
) -> List[Tuple[str, str]]:
    """
    Create the partially spliced transcripts for one sample.

    Returns:
        List of (transcript name, sequence), in a deterministic order.
    """
    partially_spliced_transcripts = []
    for transcript_id_stable in transcript_ids_stable:
        reference_transcript = reference_transcripts[transcript_id_stable]
        introns = get_introns(reference_transcript['exons'])
        rng = random.Random('%i:%s:%s' % (RANDOM_SEED, sample_id, transcript_id_stable))
        retained_intron_indices_list = sample_retained_intron_indices(
            num_introns=len(introns),
            num_samples=NUM_PARTIALLY_SPLICED_TRANSCRIPTS,
            rng=rng
        )
        assert len(retained_intron_indices_list) == NUM_PARTIALLY_SPLICED_TRANSCRIPTS, \
            '%s has %i intron(s), too few for %i distinct partially spliced transcripts' % \
            (transcript_id_stable, len(introns), NUM_PARTIALLY_SPLICED_TRANSCRIPTS)

        for i, retained_intron_indices in enumerate(retained_intron_indices_list):
            sequence = get_transcript_sequence(
                fasta=fasta,
                reference_transcript=reference_transcript,
                retained_intron_indices=retained_intron_indices
            )
            name = '%s-%s-partially-spliced-%i' % (sample_id, transcript_id_stable, i + 1)
            partially_spliced_transcripts.append((name, sequence))
            print('\t%s: %s bp, retained introns %s' % (
                name,
                len(sequence),
                ', '.join('%s:%i-%i' % (
                    reference_transcript['chromosome'], introns[j][0], introns[j][1]
                ) for j in retained_intron_indices)
            ))
    return partially_spliced_transcripts


if __name__ == '__main__':
    os.makedirs(TRANSCRIPT_DIR, exist_ok=True)

    # Step 1. Load genome data
    fasta = pysam.FastaFile(GENOME_FASTA_FILE)

    # Step 2. Load GENCODE
    gencode = Gencode(
        gtf_file=GENCODE_GTF_FILE,
        version='v45',
        species='human',
        levels=[1, 2],
        types=['protein_coding']
    )

    # Step 3. Build the reference transcript dictionary
    reference_transcripts = build_reference_transcript_dict(
        gencode=gencode,
        transcript_ids_stable=sorted({
            transcript_id_stable
            for transcript_ids_stable in REFERENCE_TRANSCRIPT_IDS.values()
            for transcript_id_stable in transcript_ids_stable
        })
    )

    # Step 4. Create the transcript files
    for sample_id in sorted(REFERENCE_TRANSCRIPT_IDS.keys()):
        print('%s:' % sample_id)
        data = {
            'transcript_id': [],
            'num_sense': [],
            'num_antisense': [],
            'sequence': []
        }

        # The mature transcripts come first. pbsim3 indexes its output by row order
        # (<sample>_chunk_%04d), so the partially spliced transcripts are appended rather than
        # interleaved, leaving the read names of the existing simulated data unchanged.
        with pysam.FastxFile('%s/%s.fasta' % (FASTA_DIR, sample_id)) as fh:
            for record in fh:
                data['transcript_id'].append(record.name)
                data['num_sense'].append(NUM_SENSE_READS)
                data['num_antisense'].append(NUM_ANTISENSE_READS)
                data['sequence'].append(record.sequence)
        num_mature = len(data['transcript_id'])

        for name, sequence in create_partially_spliced_transcripts(
                fasta=fasta,
                reference_transcripts=reference_transcripts,
                sample_id=sample_id,
                transcript_ids_stable=REFERENCE_TRANSCRIPT_IDS[sample_id]
        ):
            data['transcript_id'].append(name)
            data['num_sense'].append(NUM_SENSE_READS_PARTIALLY_SPLICED)
            data['num_antisense'].append(NUM_ANTISENSE_READS)
            data['sequence'].append(sequence)

        df = pd.DataFrame(data)
        df.to_csv(
            '%s/%s.transcript' % (TRANSCRIPT_DIR, sample_id),
            sep='\t',
            index=False,
            header=False
        )
        print('%s: %d transcripts (%d mature, %d partially spliced)' % (
            sample_id, len(df), num_mature, len(df) - num_mature
        ))
