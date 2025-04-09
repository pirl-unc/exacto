import argparse
import pandas as pd
from typing import List


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a BED file of regions to exclude for Exacto peptide variant calling."
    )

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "-i",
        dest="gtf_file",
        type=str,
        required=True,
        help="GENCODE GTF file."
    )
    parser_required.add_argument(
        "-g",
        dest="gene_list_file",
        type=str,
        required=True,
        help="IMGT gene list file."
    )
    parser_required.add_argument(
        "-s",
        dest="species",
        type=str,
        required=True,
        help="IMGT gene list species."
    )
    parser_required.add_argument(
        "-o",
        dest="output_bed_file",
        type=str,
        required=True,
        help="Output BED file."
    )
    return parser.parse_args()


def read_imgt_gene_list(gene_list_file: str, species: str) -> List[str]:
    gene_names = set()
    with open(gene_list_file, "r") as f:
        for line in f:
            elements = line.strip().split(";")
            curr_species = elements[0].lower()
            curr_gene_name = elements[1].lower()
            if curr_species == species:
                gene_names.add(curr_gene_name)
    return list(gene_names)


def run():
    # Step 1. Read the input arguments
    args = parse_args()
    args.species = args.species.lower()

    # Step 2. Read the IMGT gene list
    gene_names = read_imgt_gene_list(args.gene_list_file, args.species)
    print('%i gene names in IMGT for %s' % (len(gene_names), args.species))

    # Step 3. Read the GENCODE GTF file
    data = {
        'chromosome': [],
        'start': [],
        'end': [],
        'gene_name': []
    }
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    df_gencode_ = pd.read_csv(args.gtf_file, sep="\t", comment="#", names=columns)
    for index,row in df_gencode_[df_gencode_['feature'] == 'gene'].iterrows():
        chromosome = str(row['seqname'])
        start = int(row['start'])
        end = int(row['end'])
        attributes = row['attribute'].split(';')
        for attribute in attributes:
            attribute = attribute.lstrip()
            attribute_elements = attribute.split(" ")
            if len(attribute_elements) == 2:
                attribute_key = attribute_elements[0]
                attribute_value = attribute_elements[1]
                attribute_value = attribute_value.replace("(I)","1")
                attribute_value = attribute_value.replace("(II)","2")
                attribute_value = attribute_value.replace("(III)","3")
                attribute_value = attribute_value.replace("(IV)", "4")
                attribute_value = attribute_value.replace("(V)", "5")
                attribute_value = attribute_value.replace("(VI)", "6")
                attribute_value = attribute_value.replace("(VII)", "7")
                attribute_value = attribute_value.replace("(VIII)", "8")
                attribute_value = attribute_value.replace("(IX)", "9")
                attribute_value = attribute_value.replace("(X)", "10")
                attribute_value = attribute_value.lower()
                if attribute_key == 'gene_name':
                    data['chromosome'].append(chromosome)
                    data['start'].append(start)
                    data['end'].append(end)
                    data['gene_name'].append(attribute_value.replace('"',''))
    df_gencode = pd.DataFrame(data)
    df_gencode.drop_duplicates(inplace=True)

    # Step 4. Get the IMGT gene coordinates
    df_output = df_gencode[df_gencode['gene_name'].isin(gene_names)]
    print('%i IMGT gene names identified in GENCODE for %s' % (len(df_output['gene_name'].unique()), args.species))

    # Step 5. Write to TSV file
    df_output.to_csv(args.output_bed_file, sep="\t", index=False, header=False)
