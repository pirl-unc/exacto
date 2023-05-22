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


"""
The purpose of this python3 script is to implement the GraphGenome dataclasses.
"""


import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import re
from dataclasses import dataclass
from typing import List
from .constants import TranslocationOrientations, VariantTypes
from .nucleotide_sequence import NucleotideSequence



@dataclass(frozen=True)
class ReferenceNode:
    chromosome: str
    start: int
    end: int
    sequence: NucleotideSequence


@dataclass(frozen=True)
class VariantNode:
    source_id: str
    sample_id: str
    variant_call_id: str
    phase_block_id: str
    clone_id: str
    variant_type: str
    sequence: NucleotideSequence    # alternate allele for SnvNode
                                    # insertion sequence for InsertionNode
                                    # empty sequence for DeletionNode
                                    # empty sequence for DuplicationNode
                                    # empty sequence for InversionNode
                                    # empty sequence for TranslocationNode


@dataclass(frozen=True)
class DeletionNode(VariantNode):
    chromosome: str
    start: int
    end: int


@dataclass(frozen=True)
class DuplicationNode(VariantNode):
    chromosome: str
    start: int
    end: int


@dataclass(frozen=True)
class InsertionNode(VariantNode):
    chromosome: str
    position: int                   # insertion sequence is inserted immediately after this position


@dataclass(frozen=True)
class InversionNode(VariantNode):
    chromosome: str
    start: int
    end: int


@dataclass(frozen=True)
class SnvNode(VariantNode):
    chromosome: str
    position: int
    reference_allele: str
    alternate_allele: str


@dataclass(frozen=True)
class TranslocationNode(VariantNode):
    chromosome_1: str
    position_1: str
    chromosome_2: str
    position_2: str
    alternate_allele: str           # e.g. 'N[chr1:1000000['
    orientation: str                # e.g. 't[p['
    t_chromosome: str
    t_position: int
    p_chromosome: str
    p_position: int


@dataclass
class GrapheGenome:
    graph: nx.MultiDiGraph
    df_reference_nodes: pd.DataFrame

    def __init__(self, reference_nodes: List[ReferenceNode]):
        """
        Initializes GraphGenome object with a list of ReferenceNode objects.

        Parameters
        ----------
        reference_nodes     :   List of ReferenceNode objects.
        """
        self.graph = nx.MultiDiGraph()
        self.__node_id_counter = 0
        reference_nodes_data = {
            'node_id': [],
            'chromosome': [],
            'start': [],
            'end': []
        }
        for reference_node in reference_nodes:
            node_id = self.__get_next_node_id()
            self.graph.add_node(node_id, type=VariantTypes.REFERENCE, reference=reference_node)
            reference_nodes_data['node_id'].append(node_id)
            reference_nodes_data['chromosome'].append(reference_node.chromosome)
            reference_nodes_data['start'].append(reference_node.start)
            reference_nodes_data['end'].append(reference_node.end)
        self.df_reference_nodes = pd.DataFrame(reference_nodes_data)

    def add_deletion(
            self,
            chromosome: str,
            start: int,
            end: int,
            source_id: str = None,
            sample_id: str = None,
            variant_call_id: str = None,
            phase_block_id: str = None,
            clone_id: str = None
    ) -> int:
        """
        Add a deletion.

        Parameters
        ----------
        chromosome          :   Chromosome.
        start               :   Start position of deletion.
        end                 :   End position of deletion.
        source_id           :   Source ID (e.g. patient ID).
        sample_id           :   Sample ID (e.g. tumor sample ID).
        variant_call_id     :   Variant call ID.
        phase_block_id      :   Phase block ID.
        clone_id            :   Clone ID.

        Returns
        -------
        node_id             :   Deletion node ID.
        """
        # Step 1. Find the reference node that overlaps the (start - 1) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start - 1, end=start - 1)

        # Step 2. Split the reference node at (start - 1) position
        self.split_reference_node(node_id=reference_node_id, position=start - 1)

        # Step 3. Find the reference node that overlaps the end position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end, end=end)

        # Step 4. Split the reference node at end position
        self.split_reference_node(node_id=reference_node_id, position=end)

        # Step 5. Add deletion node
        deletion_node = DeletionNode(
            chromosome=chromosome,
            start=start,
            end=end,
            source_id=source_id,
            sample_id=sample_id,
            variant_call_id=variant_call_id,
            phase_block_id=phase_block_id,
            clone_id=clone_id,
            variant_type=VariantTypes.DELETION,
            sequence=NucleotideSequence(sequence='')
        )
        del_node_id = self.__get_next_node_id()
        self.graph.add_node(del_node_id, type=VariantTypes.DELETION, variant=deletion_node)

        # Step 6. Add an edge from the reference node at (start - 1) position to the deletion node
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start - 1, end=start - 1)
        self.graph.add_edge(reference_node_id, del_node_id, marker=False)

        # Step 7. Add an edge from the deletion node to the reference node at end + 1
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end + 1, end=end + 1)
        self.graph.add_edge(del_node_id, reference_node_id, marker=False)
        return del_node_id

    def add_duplication(
            self,
            chromosome: str,
            start: int,
            end: int,
            source_id: str = None,
            sample_id: str = None,
            variant_call_id: str = None,
            phase_block_id: str = None,
            clone_id: str = None
    ) -> int:
        """
        Add a duplication.

        Parameters
        ----------
        chromosome          :   Chromosome.
        start               :   Start position of duplication.
        end                 :   End position of duplication.
        source_id           :   Source ID (e.g. patient ID).
        sample_id           :   Sample ID (e.g. tumor sample ID).
        variant_call_id     :   Variant call ID.
        phase_block_id      :   Phase block ID.
        clone_id            :   Clone ID.

        Returns
        -------
        node_id             :   Duplication node ID.
        """
        # Step 1. Find and split the reference node that overlaps the (start - 1) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start - 1, end=start - 1)
        self.split_reference_node(node_id=reference_node_id, position=start - 1)

        # Step 2. Find and split the reference node that overlaps the (end) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end, end=end)
        self.split_reference_node(node_id=reference_node_id, position=end)

        # Step 3. Add duplication node
        duplication_node = DuplicationNode(
            chromosome=chromosome,
            start=start,
            end=end,
            source_id=source_id,
            sample_id=sample_id,
            variant_call_id=variant_call_id,
            phase_block_id=phase_block_id,
            clone_id=clone_id,
            variant_type=VariantTypes.DELETION,
            sequence=NucleotideSequence(sequence='')
        )
        dup_node_id = self.__get_next_node_id()
        self.graph.add_node(dup_node_id, type=VariantTypes.DUPLICATION, variant=duplication_node)

        # Step 4. Add an edge from the reference node at (start - 1) position to the duplication node
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start - 1, end=start - 1)
        self.graph.add_edge(reference_node_id, dup_node_id, marker=True)

        # Step 5. Add an edge from the duplication node to the reference node at (end + 1)
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end + 1, end=end + 1)
        self.graph.add_edge(dup_node_id, reference_node_id, marker=True)
        return dup_node_id

    def add_insertion(
            self,
            chromosome: str,
            position: str,
            sequence: str,
            source_id: str = None,
            sample_id: str = None,
            variant_call_id: str = None,
            phase_block_id: str = None,
            clone_id: str = None
    ) -> int:
        """
        Add an insertion.

        Parameters
        ----------
        chromosome          :   Chromosome.
        position            :   Insertion position (insertion sequence is inserted immediately after this position).
        sequence            :   Insertion sequence.
        source_id           :   Source ID (e.g. patient ID).
        sample_id           :   Sample ID (e.g. tumor sample ID).
        variant_call_id     :   Variant call ID.
        phase_block_id      :   Phase block ID.
        clone_id            :   Clone ID.

        Returns
        -------
        node_id             :   Insertion node ID.
        """
        # Step 1. Find the reference node that overlaps the INS position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position, end=position)

        # Step 2. Split the reference node at position
        self.split_reference_node(node_id=reference_node_id, position=position)

        # Step 3. Add insertion node
        insertion_node = InsertionNode(
            chromosome=chromosome,
            position=position,
            source_id=source_id,
            sample_id=sample_id,
            variant_call_id=variant_call_id,
            phase_block_id=phase_block_id,
            clone_id=clone_id,
            variant_type=VariantTypes.INSERTION,
            sequence=NucleotideSequence(sequence=sequence)
        )
        ins_node_id = self.__get_next_node_id()
        self.graph.add_node(ins_node_id, type=VariantTypes.INSERTION, variant=insertion_node)

        # Step 4. Add an edge from the reference node at position to the insertion node
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position, end=position)
        self.graph.add_edge(reference_node_id, ins_node_id, marker=False)

        # Step 5. Add an edge from the INS node to the reference node at (position + 1)
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position + 1, end=position + 1)
        self.graph.add_edge(ins_node_id, reference_node_id, marker=False)
        return ins_node_id

    def add_inversion(
            self,
            chromosome: str,
            start: int,
            end: int,
            source_id: str = None,
            sample_id: str = None,
            variant_call_id: str = None,
            phase_block_id: str = None,
            clone_id: str = None
    ) -> int:
        """
        Add an inversion.

        Parameters
        ----------
        chromosome          :   Chromosome.
        start               :   Start position of inversion.
        end                 :   End position of inversion.
        source_id           :   Source ID (e.g. patient ID).
        sample_id           :   Sample ID (e.g. tumor sample ID).
        variant_call_id     :   Variant call ID.
        phase_block_id      :   Phase block ID.
        clone_id            :   Clone ID.

        Returns
        -------
        node_id             :   Inversion node ID.
        """
        # Step 1. Find the reference node that overlaps the (start - 1) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start - 1, end=start - 1)

        # Step 2. Split the reference at (start - 1) position
        self.split_reference_node(node_id=reference_node_id, position=start - 1)

        # Step 3. Find the reference node that overlaps the (start) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start, end=start)

        # Step 4. Split the reference at (start) position
        self.split_reference_node(node_id=reference_node_id, position=start)

        # Step 5. Find the reference node that overlaps the (end - 1) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end - 1, end=end - 1)

        # Step 6. Split the reference at (end - 1) position
        self.split_reference_node(node_id=reference_node_id, position=end - 1)

        # Step 7. Find the reference node that overlaps the (end) position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end, end=end)

        # Step 8. Split the reference at (end) position
        self.split_reference_node(node_id=reference_node_id, position=end)

        # Step 9. Add inversion node
        inversion_node = InversionNode(
            chromosome=chromosome,
            start=start,
            end=end,
            source_id=source_id,
            sample_id=sample_id,
            variant_call_id=variant_call_id,
            phase_block_id=phase_block_id,
            clone_id=clone_id,
            variant_type=VariantTypes.INVERSION,
            sequence=NucleotideSequence(sequence='')
        )
        inv_node_id = self.__get_next_node_id()
        self.graph.add_node(inv_node_id, type=VariantTypes.INVERSION, variant=inversion_node)

        # Step 10. Add an edge from the reference node at (start) to the inversion node
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=start, end=start)
        self.graph.add_edge(reference_node_id, inv_node_id, marker=True)

        # Step 11. Add an edge from the inversion node to the reference node at (end)
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=end, end=end)
        self.graph.add_edge(inv_node_id, reference_node_id, marker=True)
        return inv_node_id

    def add_snv(
            self,
            chromosome: str,
            position: int,
            reference_allele: str,
            alternate_allele: str,
            source_id: str = None,
            sample_id: str = None,
            variant_call_id: str = None,
            phase_block_id: str = None,
            clone_id: str = None
    ) -> int:
        """
        Add a single-nucleotide variant.

        Parameters
        ----------
        chromosome          :   Chromosome.
        position            :   SNV position.
        reference_allele    :   Reference allele.
        alternate_allele    :   Alternate allele.
        source_id           :   Source ID (e.g. patient ID).
        sample_id           :   Sample ID (e.g. tumor sample ID).
        variant_call_id     :   Variant call ID.
        phase_block_id      :   Phase block ID.
        clone_id            :   Clone ID.

        Returns
        -------
        node_id             :   SNV node ID.
        """
        # Step 1. Find the reference node that overlaps the SNV position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position, end=position)

        # Step 2. Split the reference node at (position - 1)
        self.split_reference_node(node_id=reference_node_id, position=position - 1)

        # Step 3. Split the reference node at position
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position, end=position)
        self.split_reference_node(node_id=reference_node_id, position=position)

        # Step 4. Add SNV node
        snv_node = SnvNode(
            chromosome=chromosome,
            position=position,
            reference_allele=reference_allele,
            alternate_allele=alternate_allele,
            source_id=source_id,
            sample_id=sample_id,
            variant_call_id=variant_call_id,
            phase_block_id=phase_block_id,
            clone_id=clone_id,
            variant_type=VariantTypes.SINGLE_NUCLEOTIDE_VARIANT,
            sequence=alternate_allele
        )
        snv_node_id = self.__get_next_node_id()
        self.graph.add_node(snv_node_id, type=VariantTypes.SINGLE_NUCLEOTIDE_VARIANT, variant=snv_node)

        # Step 5. Add an edge from the reference node at (position - 1) to the SNV node
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position - 1, end=position - 1)
        self.graph.add_edge(reference_node_id, snv_node_id, marker=False)

        # Step 6. Add an edge from the SNV node to the reference node at (position + 1)
        reference_node_id = self.find_reference_node(chromosome=chromosome, start=position + 1, end=position + 1)
        self.graph.add_edge(snv_node_id, reference_node_id, marker=False)
        return snv_node_id

    def add_translocation(
            self,
            chromosome_1: str,
            position_1: int,
            chromosome_2: str,
            position_2: int,
            alternate_allele: str,
            source_id: str = None,
            sample_id: str = None,
            variant_call_id: str = None,
            phase_block_id: str = None,
            clone_id: str = None
    ) -> int:
        """
        Add a translocation.

        Parameters
        ----------
        chromosome_1        :   Chromosome 1.
        position_1          :   Position 1.
        chromosome_2        :   Chromosome 2.
        position_2          :   Position 2.
        alternate_allele    :   Alternate allele (e.g. 'N[chr1:1000000[')
        source_id           :   Source ID (e.g. patient ID).
        sample_id           :   Sample ID (e.g. tumor sample ID).
        variant_call_id     :   Variant call ID.
        phase_block_id      :   Phase block ID.
        clone_id            :   Clone ID.

        Returns
        -------
        node_id             :   Translocation node ID.
        """
        # Step 1. Parse alternate_allele information for translocation directionality
        if re.search("^.*\[.*\[$", alternate_allele):                  # t[p[ piece extending to the right of p is joined after t
            orientation = TranslocationOrientations.ORIENTATION_1
            alternate_allele_elements = alternate_allele.split('[')
            t = alternate_allele_elements[0]
            p = alternate_allele_elements[1]
        elif re.search("^.*\].*\]$", alternate_allele):                # t]p] reverse comp piece extending left of p is joined after t
            orientation = TranslocationOrientations.ORIENTATION_2
            alternate_allele_elements = alternate_allele.split(']')
            t = alternate_allele_elements[0]
            p = alternate_allele_elements[1]
        elif re.search("^\].*\].*$", alternate_allele):                # ]p]t piece extending to the left of p is joined before t
            orientation = TranslocationOrientations.ORIENTATION_3
            alternate_allele_elements = alternate_allele.split(']')
            t = alternate_allele_elements[2]
            p = alternate_allele_elements[1]
        elif re.search("^\[.*\[.*$", alternate_allele):                # [p[t  reverse comp piece extending right of p is joined before t
            orientation = TranslocationOrientations.ORIENTATION_4
            alternate_allele_elements = alternate_allele.split('[')
            t = alternate_allele_elements[2]
            p = alternate_allele_elements[1]
        else:
            raise Exception('Unknown ALT format to infer translocation orientation type: %s' % alternate_allele)
        if p == '%s:%i' % (chromosome_1, position_1):
            p_chromosome = chromosome_1
            p_position = position_1
            t_chromosome = chromosome_2
            t_position = position_2
        elif p == '%s:%i' % (chromosome_2, position_2):
            p_chromosome = chromosome_2
            p_position = position_2
            t_chromosome = chromosome_1
            t_position = position_1
        else:
            raise Exception('Positions for p and t could not be inferred from self.alternate_allele: %s' % alternate_allele)

        # Step 2. Depending on the translocation directionality,
        # find and split the reference node that overlaps t_position and p_position.
        # Here the goal is to split reference nodes such that the breakpoint in positions
        # t and p are single nodes encoding single nucleotide.
        if orientation == TranslocationOrientations.ORIENTATION_1:
            # Split at (t_position - 1)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position - 1, end=t_position - 1)
            self.split_reference_node(node_id=reference_node_id, position=t_position)

            # Split at (t_position)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.split_reference_node(node_id=reference_node_id, position=t_position)

            # Split at (p_position - 1)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position - 1, end=p_position - 1)
            self.split_reference_node(node_id=reference_node_id, position=p_position - 1)

            # Split at (p_position)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.split_reference_node(node_id=reference_node_id, position=p_position)
        elif orientation == TranslocationOrientations.ORIENTATION_2:
            # Split at (t_position - 1)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position - 1, end=t_position - 1)
            self.split_reference_node(node_id=reference_node_id, position=t_position - 1)

            # Split at (t_position)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.split_reference_node(node_id=reference_node_id, position=t_position)

            # Split at (p_position - 1)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position - 1, end=p_position - 1)
            self.split_reference_node(node_id=reference_node_id, position=p_position - 1)

            # Split at (p_position)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.split_reference_node(node_id=reference_node_id, position=p_position)
        elif orientation == TranslocationOrientations.ORIENTATION_3:
            # Split at (t_position)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.split_reference_node(node_id=reference_node_id, position=t_position)

            # Split at (t_position + 1)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position + 1, end=t_position + 1)
            self.split_reference_node(node_id=reference_node_id, position=t_position + 1)

            # Split at (p_position - 1)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position - 1, end=p_position - 1)
            self.split_reference_node(node_id=reference_node_id, position=p_position - 1)

            # Split at (p_position)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.split_reference_node(node_id=reference_node_id, position=p_position)
        elif orientation == TranslocationOrientations.ORIENTATION_4:
            # Split at (t_position)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.split_reference_node(node_id=reference_node_id, position=t_position)

            # Split at (t_position + 1)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position + 1, end=t_position + 1)
            self.split_reference_node(node_id=reference_node_id, position=t_position + 1)

            # Split at (p_position - 1)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position - 1, end=p_position - 1)
            self.split_reference_node(node_id=reference_node_id, position=p_position - 1)

            # Split at (p_position)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.split_reference_node(node_id=reference_node_id, position=p_position)
        else:
            raise Exception('Unknown translocation orientation: %s' % orientation)

        # Step 3. Add translocation node
        translocation_node = TranslocationNode(
            chromosome_1=chromosome_1,
            position_1=position_1,
            chromosome_2=chromosome_2,
            position_2=position_2,
            alternate_allele=alternate_allele,
            orientation=orientation,
            t_chromosome=t_chromosome,
            t_position=t_position,
            p_chromosome=p_chromosome,
            p_position=p_position,
            source_id=source_id,
            sample_id=sample_id,
            variant_call_id=variant_call_id,
            phase_block_id=phase_block_id,
            clone_id=clone_id,
            variant_type=VariantTypes.TRANSLOCATION,
            sequence=NucleotideSequence(sequence='')
        )
        tra_node_id = self.__get_next_node_id()
        self.graph.add_node(tra_node_id, type=VariantTypes.TRANSLOCATION, variant=translocation_node)

        # Step 4. Depending on the translocation directionality,
        # add edges between the translocation node and the reference nodes.
        if orientation == TranslocationOrientations.ORIENTATION_1:
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.graph.add_edge(reference_node_id, tra_node_id, marker=True)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.graph.add_edge(tra_node_id, reference_node_id, marker=True)
        elif orientation == TranslocationOrientations.ORIENTATION_2:
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.graph.add_edge(reference_node_id, tra_node_id, marker=True)
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.graph.add_edge(tra_node_id, reference_node_id, marker=True)
        elif orientation == TranslocationOrientations.ORIENTATION_3:
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.graph.add_edge(reference_node_id, tra_node_id, marker=True)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.graph.add_edge(tra_node_id, reference_node_id, marker=True)
        elif orientation == TranslocationOrientations.ORIENTATION_4:
            reference_node_id = self.find_reference_node(chromosome=p_chromosome, start=p_position, end=p_position)
            self.graph.add_edge(reference_node_id, tra_node_id, marker=True)
            reference_node_id = self.find_reference_node(chromosome=t_chromosome, start=t_position, end=t_position)
            self.graph.add_edge(tra_node_id, reference_node_id, marker=True)
        else:
            raise Exception('Unknown translocation orientation: %s' % orientation)
        return tra_node_id

    def find_all_simple_paths(self, source_node_id: str, target_node_id: str):
        """
        Finds all simple paths from source node ID to target node ID.

        Parameters
        ----------
        source_node_id  :   str
        target_node_id  :   str

        Returns
        -------
        paths           :   List of list of tuples (e.g. [[(1,2), (2,3)], [(1,3)]])
        """
        paths = nx.all_simple_paths(self.graph, source=source_node_id, target=target_node_id)
        paths_list = []
        for path in map(nx.utils.pairwise, paths):
            paths_list.append(list(path))
        return paths_list

    def find_reference_node(self, chromosome, start, end):
        """
        Finds a reference node that matches the query position.

        Parameters
        ----------
        chromosome  :   Chromosome.
        start       :   Start.
        end         :   End.

        Returns
        -------
        node_id     :   Reference node ID.
        """
        df_matched = self.df_reference_nodes[
            (self.df_reference_nodes['chromosome'] == chromosome) &
            (self.df_reference_nodes['start'] <= end) &
            (self.df_reference_nodes['end'] >= start)
        ]
        if len(df_matched) == 1:
            return df_matched['node_id'].values.tolist()[0]
        elif len(df_matched) > 1:
            raise Exception('More than 1 reference node was found for %s:%i-%i.' % (chromosome, start, end))
        else:
            raise Exception('No reference node could be found for %s:%i-%i.' % (chromosome, start, end))

    def _minimizers(self, window, k):
        """


        Parameters
        ----------
        window
        k

        Returns
        -------

        """

    def plot(
            self,
            reference_hex: str = '#DCDCDC',
            snv_hex: str = '#C8A200',
            deletion_hex: str = '#4667EA',
            insertion_hex: str = '#EF5358',
            duplication_hex: str = '#ff9933',
            inversion_hex: str = '#478E0E',
            translocation_hex: str = '#AD52BF',
            width: float = 16,
            height: float = 9,
            font_size: float = 12,
            node_size: float = 600,
            arroww_size: float = 14,
            line_widths: float = 2
    ) -> plt.figure:
        # Step 1. Assign node colors
        node_colors = []
        for node_id, node_data in self.graph.nodes(data=True):
            if node_data['type'] == VariantTypes.REFERENCE:
                node_colors.append(reference_hex)
            elif node_data['type'] == VariantTypes.SINGLE_NUCLEOTIDE_VARIANT:
                node_colors.append(snv_hex)
            elif node_data['type'] == VariantTypes.DELETION:
                node_colors.append(deletion_hex)
            elif node_data['type'] == VariantTypes.INSERTION:
                node_colors.append(insertion_hex)
            elif node_data['type'] == VariantTypes.DUPLICATION:
                node_colors.append(duplication_hex)
            elif node_data['type'] == VariantTypes.INVERSION:
                node_colors.append(inversion_hex)
            elif node_data['type'] == VariantTypes.TRANSLOCATION:
                node_colors.append(translocation_hex)

        # Step 2. Assign edge styles
        markers = [self.graph[u][v][0]['marker'] for u, v in self.graph.edges()]
        edge_styles = []
        for marker in markers:
            if marker:
                edge_styles.append('dashed')
            else:
                edge_styles.append('solid')

        # Step 2. Draw graph
        fig, ax = plt.subplots(figsize=(width, height))
        nx.draw_networkx(
            self.graph,
            ax=ax,
            node_color=node_colors,
            style=edge_styles,
            font_size=font_size,
            node_size=node_size,
            arrowsize=arroww_size,
            linewidths=line_widths,
            pos=nx.nx_pydot.graphviz_layout(self.graph)
        )

        # Step 3. Add a legend
        patch_ref = mpatches.Patch(color=reference_hex, label=VariantTypes.FULL_NAMES[VariantTypes.REFERENCE])
        patch_snv = mpatches.Patch(color=snv_hex, label=VariantTypes.FULL_NAMES[VariantTypes.SINGLE_NUCLEOTIDE_VARIANT])
        patch_del = mpatches.Patch(color=deletion_hex, label=VariantTypes.FULL_NAMES[VariantTypes.DELETION])
        patch_ins = mpatches.Patch(color=insertion_hex, label=VariantTypes.FULL_NAMES[VariantTypes.INSERTION])
        patch_dup = mpatches.Patch(color=duplication_hex, label=VariantTypes.FULL_NAMES[VariantTypes.DUPLICATION])
        patch_inv = mpatches.Patch(color=inversion_hex, label=VariantTypes.FULL_NAMES[VariantTypes.INVERSION])
        patch_tra = mpatches.Patch(color=translocation_hex, label=VariantTypes.FULL_NAMES[VariantTypes.TRANSLOCATION])
        plt.legend(
            handles=[patch_ref, patch_snv, patch_del, patch_ins, patch_dup, patch_inv, patch_tra],
            fontsize=font_size
        )
        return fig

    def split_reference_node(self, node_id, position):
        """
        Splits a node into two: [:position] and [position+1:]

        Parameters
        ----------
        node_id     :   Reference node ID.
        position    :   Position to split.

        Returns
        -------

        """
        # Step 1. Fetch the reference node object
        reference_node = self.graph.nodes[node_id]['reference']

        # Step 2. Find all incoming nodes to the reference node
        incoming_node_ids = list(self.graph.predecessors(node_id))
        incoming_node_edge_makers = []
        for incoming_node_id in incoming_node_ids:
            incoming_node_edge_makers.append(self.graph.get_edge_data(incoming_node_id, node_id)[0]['marker'])

        # Step 3. Find al outgoing nodes from the reference node
        outgoing_node_ids = list(self.graph.successors(node_id))
        outgoing_node_edge_makers = []
        for outgoing_node_id in outgoing_node_ids:
            outgoing_node_edge_makers.append(self.graph.get_edge_data(node_id, outgoing_node_id)[0]['marker'])

        # Step 4. Create two reference node objects (to split them)
        reference_node_1 = ReferenceNode(
            chromosome=reference_node.chromosome,
            start=reference_node.start,
            end=position,
            sequence=NucleotideSequence(
                sequence=str(reference_node.sequence)[:(position - reference_node.start + 1)]
            )
        )
        reference_node_2 = ReferenceNode(
            chromosome=reference_node.chromosome,
            start=position + 1,
            end=reference_node.end,
            sequence=NucleotideSequence(
                sequence=str(reference_node.sequence)[(position - reference_node.start + 1):]
            )
        )

        # Step 5. Add the two new nodes
        reference_node_1_id = self.__get_next_node_id()
        reference_node_2_id = self.__get_next_node_id()
        self.graph.add_node(reference_node_1_id, type=VariantTypes.REFERENCE, reference=reference_node_1)
        self.graph.add_node(reference_node_2_id, type=VariantTypes.REFERENCE, reference=reference_node_2)
        self.df_reference_nodes = self.df_reference_nodes[self.df_reference_nodes['node_id'] != node_id]
        new_reference_nodes_data = {
            'node_id': [reference_node_1_id, reference_node_2_id],
            'chromosome': [reference_node_1.chromosome, reference_node_2.chromosome],
            'start': [reference_node_1.start, reference_node_2.start],
            'end': [reference_node_1.end, reference_node_2.end]
        }
        self.df_reference_nodes = pd.concat([self.df_reference_nodes, pd.DataFrame(new_reference_nodes_data)])

        # Step 6. Add edges from all incoming nodes to reference_node to reference_node_1
        for i in range(0, len(incoming_node_ids)):
            incoming_node_id = incoming_node_ids[i]
            self.graph.add_edge(incoming_node_id, reference_node_1_id, marker=incoming_node_edge_makers[i])

        # Step 7. Add an edge from reference_node_1 to reference_node_2
        self.graph.add_edge(reference_node_1_id, reference_node_2_id, marker=False)

        # Step 8. Add edges from reference_node_2 to all outgoing nodes from reference_node
        for i in range(0, len(outgoing_node_ids)):
            outgoing_node_id = outgoing_node_ids[i]
            self.graph.add_edge(reference_node_2_id, outgoing_node_id, marker=outgoing_node_edge_makers[i])

        # Step 9. Remove original reference node
        self.graph.remove_node(node_id)

    def __get_next_node_id(self) -> int:
        self.__node_id_counter += 1
        return self.__node_id_counter
