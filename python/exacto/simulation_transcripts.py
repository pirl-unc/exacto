# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# import pysam
# import pandas as pd
# from typing import List
# from exacto.gene import Gene
# from exacto.transcript import Transcript
# from exacto.exon import Exon
#
#
# def simulate_transcripts(
#         df_ref_genes: pd.DataFrame,
#         df_ref_transcripts: pd.DataFrame,
#         df_ref_exons: pd.DataFrame,
#         ref_fasta: pysam.FastaFile,
#         df_germline_variants: pd.DataFrame,
#         df_somatic_variants: pd.DataFrame
#     ) -> List[Gene]:
#     """
#     Simulates transcript sequences.
#
#     Args
#     ----
#     df_ref_genes        :   DataFrame with the following columns:
#                             'gene_id'
#                             'gene_name'
#                             'gene_type'
#                             'gene_chrom'
#                             'gene_start'
#                             'gene_end'
#                             'gene_strand'
#     df_ref_transcripts  :   DataFrame with the following columns:
#                             'gene_id'
#                             'transcript_id'
#                             'transcript_name'
#                             'transcript_type'
#                             'transcript_chrom'
#                             'transcript_start'
#                             'transcript_end'
#                             'transcript_strand'
#     df_ref_exons        :   DataFrame with the following columns:
#                             'gene_id'
#                             'transcript_id'
#                             'exon_id'
#                             'exon_number'
#                             'exon_chrom'
#                             'exon_start'
#                             'exon_end'
#     ref_fasta           :   pysam.FastaFile object.
#
#     Returns
#     -------
#     """
#     # Step 1. Load reference genes, transcripts, and exons
#     genes = [] # List of Gene class objects
#     for _, curr_gene_row in df_ref_genes.iterrows():
#         curr_gene_id = curr_gene_row['gene_id']
#         curr_gene_chromosome = curr_gene_row['gene_chrom']
#         curr_gene_start = curr_gene_row['gene_start']
#         curr_gene_end = curr_gene_row['gene_end']
#         curr_gene_strand = curr_gene_row['gene_strand']
#         gene = Gene(
#             ref_id=curr_gene_id,
#             ref_chromosome=curr_gene_chromosome,
#             ref_start=curr_gene_start,
#             ref_end=curr_gene_end,
#             ref_strand=curr_gene_strand
#         )
#         df_curr_transcripts = df_ref_transcripts.loc[df_ref_transcripts['gene_id'] == curr_gene_id,:]
#         for _, curr_transcript_row in df_curr_transcripts.iterrows():
#             curr_transcript_id = curr_transcript_row['transcript_id']
#             curr_transcript_chromosome = curr_transcript_row['transcript_chrom']
#             curr_transcript_start = curr_transcript_row['transcript_start']
#             curr_transcript_end = curr_transcript_row['transcript_end']
#             curr_transcript_type = curr_transcript_row['transcript_type']
#             curr_transcript_strand = curr_transcript_row['transcript_strand']
#             transcript = Transcript(
#                 ref_id=curr_transcript_id,
#                 ref_chromosome=curr_transcript_chromosome,
#                 ref_start=curr_transcript_start,
#                 ref_end=curr_transcript_end,
#                 ref_type=curr_transcript_type,
#                 ref_strand=curr_transcript_strand
#             )
#             df_curr_exons = df_ref_exons.loc[df_ref_exons['transcript_id'] == curr_transcript_id,:]
#             df_curr_exons = df_curr_exons.sort_values(['exon_number'], ascending=True)
#             for _, curr_exon in df_curr_exons.iterrows():
#                 curr_exon_id = curr_exon['exon_id']
#                 curr_exon_start = curr_exon['exon_start']
#                 curr_exon_end = curr_exon['exon_end']
#                 curr_exon_sequence = ref_fasta.fetch(curr_transcript_chromosome,
#                                                      curr_exon_start - 1,
#                                                      curr_exon_end)
#                 exon = Exon(
#                     ref_id=curr_exon_id,
#                     ref_start=curr_exon_start,
#                     ref_end=curr_exon_end,
#                     ref_sequence=curr_exon_sequence
#                 )
#                 transcript.insert_exon(exon=exon)
#             gene.add_transcript(transcript=transcript)
#         genes.append(gene)
#
#     # Step 2. Apply germline variants to each reference transcript
#
#     # Step 3. Apply somatic variants to each germline transcript
