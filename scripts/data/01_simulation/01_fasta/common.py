import random
from dataclasses import dataclass, replace
from typing import List

import pandas as pd


def create_fasta_file(
        sequences: List[str],
        sequence_names: List[str],
        output_fasta_file: str
):
    """
    Create a FASTA file.

    Parameters:
        sequences               :   List of sequences.
        output_fasta_file       :   Output FASTA file.
    """
    with open(output_fasta_file, 'w') as f:
        for i, sequence in enumerate(sequences):
            f.write('>%s\n' % sequence_names[i])
            f.write(sequence + '\n')


def generate_random_sequence(k: int) -> str:
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(nucleotides, k=k))


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


_COMPLEMENT = str.maketrans(
    'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
    'TGCAYRSWMKVHDBNtgcayrswmkvhdbn'
)


_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def _translate(seq, frame):
    return ''.join(
        _CODON_TABLE.get(seq[i:i + 3], 'X')
        for i in range(frame, len(seq) - 2, 3)
    )


def longest_orf_peptide(seq, both_strands=False, require_stop=True):
    """Return the peptide of the longest ORF in a nucleotide sequence.

    An ORF is ATG -> in-frame stop. The peptide starts with 'M' and excludes
    the stop. Non-ACGT codons translate to 'X'. Returns '' if no ORF is found.

    both_strands : also scan the reverse complement (6-frame). Set False for
                   strand-resolved transcripts where only the + strand is sense.
    require_stop : only count ORFs that terminate in a stop codon; if False,
                   an ORF may run off the end of the sequence.
    """
    seq = seq.upper()
    strands = [seq, seq.translate(_COMPLEMENT)[::-1]] if both_strands else [seq]

    best = ''
    for s in strands:
        for frame in range(3):
            protein = _translate(s, frame)
            start = 0
            while start < len(protein):
                stop = protein.find('*', start)
                segment = protein[start:] if stop == -1 else protein[start:stop]
                m = segment.find('M')
                if m != -1 and (stop != -1 or not require_stop):
                    if len(segment) - m > len(best):
                        best = segment[m:]
                if stop == -1:
                    break
                start = stop + 1
    return best


# Transcript model alignments ground truth.
#
# A TranscriptStructure records how a simulated transcript is put together from the genome, so the
# alignment model of the transcript can be written next to its sequence. Its operations mirror the
# string operations of the create_rna scripts: + strand fetches of exons and segments, substituted
# and inserted bases, concatenation (+) and reverse_complement(). get_sequence() spells the
# transcript, and write_transcript_model_alignments_ground_truth() refuses a structure that does not
# spell the sequence the script wrote.
#
# The rows have the columns of the transcript model alignments table of exacto call-rna-vars
# (<prefix>_exacto_assembled_transcript_model_alignments.tsv), in transcript (read) order, for the
# transcript aligned as a read:
#   base rows      a run of matched, substituted or inserted bases. read_start..read_end is
#                  inclusive. A match spans its reference bases (I, I); a substitution and an
#                  insertion are given by the reference bases either side (D, U).
#   event rows     a join between two runs that do not continue each other: a deletion inside
#                  one exon, a splice (a forward jump on one chromosome and strand, at most
#                  MAX_INTRON_LENGTH bases), or a breakpoint (any other jump). read_start and
#                  read_end are the last base before and the first base after it.
#   positions      two sides in ascending position order. A side where its run ends going up the
#                  chromosome is D, where it starts is U.
#   sequence       + strand bases, upper case; inserted bases at a breakpoint go in its row.
#   context        a base is exonic in an exon of its source transcript, intronic elsewhere in
#                  that transcript, intergenic outside a transcript. A splice is canonical when
#                  it skips exactly an intron of a source transcript, fusion when it joins two
#                  genes, noncanonical otherwise. A breakpoint is backsplicing when bases before
#                  and after it share a position and strand, fusion when it joins two genes,
#                  noncanonical on one chromosome otherwise.
#   skipped        the exonic bases of a source transcript an event skips that the transcript holds
#                  nowhere: between the two sides of a forward join within one transcript; after
#                  the left side and before the right side of a join between two transcripts. One
#                  run per exon, chr:start:S:strand|chr:end:S:strand|gene:transcript:exon|gene:transcript:exon,
#                  in position order, separated by ';'.
#   reference_gene_name, reference_transcript_id
#                  the genes and transcripts the bases come from, in transcript order, separated
#                  by ';'.
# The genes, transcripts and exons are those the script took the bases from, not those an
# annotation lookup would find at the position.

MAX_INTRON_LENGTH = 200000  # minimap2 -G default: a longer forward jump is split into records

TRANSCRIPT_MODEL_ALIGNMENT_COLUMNS = [
    'assembled_transcript_name', 'reference_gene_name', 'reference_transcript_id', 'index',
    'read_start', 'read_end', 'sequence', 'type', 'kind', 'context',
    'chromosome_1', 'position_1', 'operation_1', 'strand_1',
    'chromosome_2', 'position_2', 'operation_2', 'strand_2',
    'reference_gene_id_1', 'reference_transcript_id_1', 'reference_exon_id_1',
    'reference_gene_id_2', 'reference_transcript_id_2', 'reference_exon_id_2',
    'skipped'
]


@dataclass(frozen=True)
class _Segment:
    chromosome: str
    start: int              # 1-based, inclusive
    end: int                # 1-based, inclusive
    strand: str             # the strand the transcript reads: '+' or '-'
    gene_id: str
    transcript_id: str
    exon_id: str
    bases: str              # + strand bases: the reference, or the bases put in its place
    substituted: bool


@dataclass(frozen=True)
class _Insertion:
    bases: str              # + strand orientation of the inserted bases
    strand: str


class TranscriptStructure:
    def __init__(self, pieces=()):
        self.pieces = tuple(pieces)

    @staticmethod
    def segment(fasta, chromosome, start, end, gene_id='', transcript_id='', exon_id=''):
        """The + strand bases chromosome:start-end (1-based, inclusive)."""
        start, end = int(start), int(end)
        assert start <= end, '%s:%i-%i is empty' % (chromosome, start, end)
        bases = str(fasta.fetch(chromosome, start - 1, end))
        return TranscriptStructure([_Segment(chromosome, start, end, '+', gene_id, transcript_id, exon_id, bases, False)])

    @staticmethod
    def exon(fasta, row, start=None, end=None):
        """A GENCODE exon row (df_exons), or its bases from start to end."""
        start = int(row['start']) if start is None else int(start)
        end = int(row['end']) if end is None else int(end)
        assert int(row['start']) <= start and end <= int(row['end']), '%i-%i is outside exon %s' % (start, end, row['exon_id'])
        return TranscriptStructure.segment(fasta, row['chromosome'], start, end, row['gene_id'], row['transcript_id'], row['exon_id'])

    @staticmethod
    def substitution(fasta, row, position, bases):
        """bases (+ strand) in place of the reference from position on, inside a GENCODE exon row."""
        position = int(position)
        exon = TranscriptStructure.exon(fasta, row, position, position + len(bases) - 1).pieces[0]
        assert exon.bases.upper() != bases.upper(), 'the substitution at %i changes no base' % position
        return TranscriptStructure([replace(exon, bases=bases, substituted=True)])

    @staticmethod
    def insertion(bases):
        """Bases from no reference position, in the orientation of the pieces around them."""
        return TranscriptStructure([_Insertion(bases, '+')])

    def __add__(self, other):
        return TranscriptStructure(self.pieces + other.pieces)

    def reverse_complement(self):
        flip = {'+': '-', '-': '+'}
        return TranscriptStructure(replace(piece, strand=flip[piece.strand]) for piece in reversed(self.pieces))

    def get_sequence(self):
        return ''.join(piece.bases if piece.strand == '+' else reverse_complement(piece.bases) for piece in self.pieces)

    def get_rows(self, name, gencode):
        """The rows of this transcript, named name, as described at the top of this section."""
        pieces = _merge_segments(self.pieces)
        segments = [piece for piece in pieces if isinstance(piece, _Segment)]
        transcript_ids = list(dict.fromkeys(s.transcript_id for s in segments if s.transcript_id))
        gene_names = []
        introns = set()
        exons = {}                  # transcript ID: [(chromosome, start, end, strand, gene ID, exon ID)] by start
        for transcript_id in transcript_ids:
            df_exons = gencode.df_exons[gencode.df_exons['transcript_id'] == transcript_id].sort_values(by=['start'])
            exons[transcript_id] = [(row['chromosome'], int(row['start']), int(row['end']), row['strand'], row['gene_id'], row['exon_id'])
                                    for _, row in df_exons.iterrows()]
            introns.update((exons[transcript_id][i][2] + 1, exons[transcript_id][i + 1][1] - 1) for i in range(len(exons[transcript_id]) - 1))
            gene_id = exons[transcript_id][0][4]
            gene_names.append(gencode.df_genes.loc[gencode.df_genes['gene_id'] == gene_id, 'name'].values[0])

        def skipped_bases(a, b, is_forward):
            """The skipped cell of the event from segment a to segment b."""
            ranges = []             # (transcript ID, first position, last position)
            if a.transcript_id and a.transcript_id == b.transcript_id:
                if is_forward:
                    low, high = sorted([_exit(a)[0], _entry(b)[0]])
                    ranges.append((a.transcript_id, low + 1, high - 1))
            else:
                if a.transcript_id:
                    ranges.append((a.transcript_id, a.end + 1, float('inf')) if a.strand == '+' else (a.transcript_id, 0, a.start - 1))
                if b.transcript_id:
                    ranges.append((b.transcript_id, 0, b.start - 1) if b.strand == '+' else (b.transcript_id, b.end + 1, float('inf')))
            runs = []
            for transcript_id, first, last in ranges:
                for chromosome, start, end, strand, gene_id, exon_id in exons[transcript_id]:
                    start, end = max(start, first), min(end, last)
                    # The exon bases in range, minus the bases the transcript holds.
                    pieces = [(start, end)] if start <= end else []
                    for segment in segments:
                        if segment.chromosome != chromosome or segment.strand != strand:
                            continue
                        pieces = [part for piece_start, piece_end in pieces for part in
                                  [(piece_start, min(piece_end, segment.start - 1)), (max(piece_start, segment.end + 1), piece_end)]
                                  if part[0] <= part[1]]
                    label = '%s:%s:%s' % (gene_id, transcript_id, exon_id)
                    runs.extend((chromosome, piece_start, piece_end, strand, label) for piece_start, piece_end in pieces)
            return ';'.join('%s:%i:S:%s|%s:%i:S:%s|%s|%s' % (chromosome, start, strand, chromosome, end, strand, label, label)
                            for chromosome, start, end, strand, label in sorted(set(runs), key=lambda run: (run[0], run[1], run[3])))

        rows = []
        def add_row(read_start, read_end, sequence, record_type, kind, context, sides, skipped=''):
            (s1, p1, o1), (s2, p2, o2) = sorted(sides, key=lambda side: side[1])
            rows.append({
                'assembled_transcript_name': name,
                'reference_gene_name': ';'.join(gene_names),
                'reference_transcript_id': ';'.join(transcript_ids),
                'index': len(rows),
                'read_start': read_start,
                'read_end': read_end,
                'sequence': sequence.upper(),
                'type': record_type,
                'kind': kind,
                'context': context,
                'chromosome_1': s1.chromosome, 'position_1': p1, 'operation_1': o1, 'strand_1': s1.strand,
                'chromosome_2': s2.chromosome, 'position_2': p2, 'operation_2': o2, 'strand_2': s2.strand,
                'reference_gene_id_1': s1.gene_id, 'reference_transcript_id_1': s1.transcript_id, 'reference_exon_id_1': s1.exon_id,
                'reference_gene_id_2': s2.gene_id, 'reference_transcript_id_2': s2.transcript_id, 'reference_exon_id_2': s2.exon_id,
                'skipped': skipped
            })

        read_position = 0
        previous = None             # (index in pieces, segment, its last read position)
        insertion = None            # (insertion, its first read position)
        for i, piece in enumerate(pieces):
            if isinstance(piece, _Insertion):
                assert previous is not None and insertion is None, '%s: an insertion must follow a segment' % name
                insertion = (piece, read_position)
                read_position += len(piece.bases)
                continue
            if previous is not None:
                j, a, a_last = previous
                exit_side, entry_side = (a,) + _exit(a), (piece,) + _entry(piece)
                if _continues(a, piece):
                    if insertion is not None:
                        inserted, first = insertion
                        add_row(first, read_position - 1, inserted.bases, 'base', 'insertion', _context(a), [exit_side, entry_side])
                elif _is_splice(a, piece) and a.exon_id and a.exon_id == piece.exon_id:
                    add_row(a_last, read_position, insertion[0].bases if insertion else '', 'event', 'deletion', '', [exit_side, entry_side],
                            skipped_bases(a, piece, True))
                elif _is_splice(a, piece):
                    skipped = (min(exit_side[1], entry_side[1]) + 1, max(exit_side[1], entry_side[1]) - 1)
                    context = 'canonical' if skipped in introns else ('fusion' if _joins_genes(a, piece) else 'noncanonical')
                    add_row(a_last, read_position, insertion[0].bases if insertion else '', 'event', 'splicing', context, [exit_side, entry_side],
                            skipped_bases(a, piece, True))
                else:
                    before = [s for s in pieces[:j + 1] if isinstance(s, _Segment)]
                    after = [s for s in pieces[i:] if isinstance(s, _Segment)]
                    if any(_overlaps(x, y) for x in before for y in after):
                        context = 'backsplicing'
                    elif _joins_genes(a, piece):
                        context = 'fusion'
                    elif a.chromosome == piece.chromosome:
                        context = 'noncanonical'
                    else:
                        context = ''
                    add_row(a_last, read_position, insertion[0].bases if insertion else '', 'event', 'breakpoint', context, [exit_side, entry_side],
                            skipped_bases(a, piece, False))
            read_start, read_end = read_position, read_position + piece.end - piece.start
            if piece.substituted:
                add_row(read_start, read_end, piece.bases, 'base', 'mismatch', _context(piece),
                        [(piece, piece.start - 1, 'D'), (piece, piece.end + 1, 'U')])
            else:
                add_row(read_start, read_end, piece.bases, 'base', 'match', _context(piece),
                        [(piece, piece.start, 'I'), (piece, piece.end, 'I')])
            previous = (i, piece, read_end)
            insertion = None
            read_position = read_end + 1
        assert insertion is None, '%s: an insertion must be followed by a segment' % name
        return rows


def _context(segment):
    return 'exonic' if segment.exon_id else ('intronic' if segment.transcript_id else 'intergenic')


def _exit(segment):
    """(position, operation) where the transcript leaves a segment."""
    return (segment.end, 'D') if segment.strand == '+' else (segment.start, 'U')


def _entry(segment):
    """(position, operation) where the transcript enters a segment."""
    return (segment.start, 'U') if segment.strand == '+' else (segment.end, 'D')


def _continues(a, b):
    """Whether b starts at the base after a, in the direction the transcript reads."""
    return a.chromosome == b.chromosome and a.strand == b.strand and \
        (b.start == a.end + 1 if a.strand == '+' else b.end == a.start - 1)


def _is_splice(a, b):
    """Whether the transcript jumps forward from a to b on one chromosome and strand, within the
    longest intron the aligner joins in one record."""
    if a.chromosome != b.chromosome or a.strand != b.strand:
        return False
    gap = b.start - a.end - 1 if a.strand == '+' else a.start - b.end - 1
    return 0 < gap <= MAX_INTRON_LENGTH


def _joins_genes(a, b):
    return a.gene_id != '' and b.gene_id != '' and a.gene_id != b.gene_id


def _overlaps(a, b):
    return a.chromosome == b.chromosome and a.strand == b.strand and a.start <= b.end and b.start <= a.end


def _merge_segments(pieces):
    """Joins neighbouring segments that continue each other with the same source and no
    substitution, as the aligner reads them as one run."""
    merged = []
    for piece in pieces:
        last = merged[-1] if merged else None
        if isinstance(piece, _Segment) and isinstance(last, _Segment) and _continues(last, piece) \
                and not last.substituted and not piece.substituted \
                and (last.gene_id, last.transcript_id, last.exon_id) == (piece.gene_id, piece.transcript_id, piece.exon_id):
            bases = last.bases + piece.bases if last.strand == '+' else piece.bases + last.bases
            merged[-1] = replace(last, start=min(last.start, piece.start), end=max(last.end, piece.end), bases=bases)
        else:
            merged.append(piece)
    return merged


def write_transcript_model_alignments_ground_truth(gencode, transcripts, output_tsv_file):
    """
    Write the transcript model alignments of each transcript.

    Parameters:
        gencode             :   vstolib Gencode the transcripts were built from (for introns).
        transcripts         :   List of (name, TranscriptStructure, sequence written to the FASTA).
        output_tsv_file     :   Output TSV file.
    """
    rows = []
    for name, structure, sequence in transcripts:
        assert structure.get_sequence() == sequence, '%s: the structure does not spell the transcript' % name
        rows.extend(structure.get_rows(name, gencode))
    pd.DataFrame(rows, columns=TRANSCRIPT_MODEL_ALIGNMENT_COLUMNS).to_csv(output_tsv_file, sep='\t', index=False)
