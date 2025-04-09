import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Format TranSigner assignments.tsv for Exacto RNA variant calling."
    )

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "-i",
        dest="tsv_file",
        type=str,
        required=True,
        help="TranSigner assignments.tsv file."
    )
    parser_required.add_argument(
        "-o",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "-t",
        dest="threshold",
        type=float,
        required=False,
        default=0.5,
        help="Probability threshold (default: 0.5). Only transcripts with this threshold value and above will be retained."
    )
    parser_optional.add_argument(
        "-s",
        dest="source",
        type=str,
        required=False,
        default='gencode',
        help="Reference FASTA source (default: 'gencode')."
    )
    return parser.parse_args()


def run():
    # Step 1. Read the input arguments
    args = parse_args()

    # Step 2. Read the data into a DataFrame
    data = {
        'read_id': [],
        'transcript_id': [],
        'probability': []
    }
    with open(args.tsv_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            first = True
            read_id = ''
            for column in columns:
                if first:
                    read_id = column
                    first = False
                else:
                    if args.source == 'gencode':
                        column = column[1:-1] # drop the ( and ) at the start and end
                        transcript_id = column.split('|')[0]
                        probability = float(column.split(',')[1])
                        if probability >= args.threshold:
                            data['read_id'].append(read_id)
                            data['transcript_id'].append(transcript_id)
                            data['probability'].append(probability)
                    else:
                        raise Exception('Unsupported reference FASTA source: %s' % args.source)
    df = pd.DataFrame(data)

    # Step 3. Apply probability threshold
    df = df[df['probability'] >= args.threshold]

    # Step 4. Write the DataFrame to a TSV file
    df.to_csv(args.output_tsv_file, sep='\t', index=False)
