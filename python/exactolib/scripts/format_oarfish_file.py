import argparse
import pandas as pd


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise Exception('Boolean value expected.')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Format Oarfish prob file for Exacto RNA variant calling."
    )

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "-i",
        dest="prob_file",
        type=str,
        required=True,
        help="Oarfish prob file."
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
        choices=['gencode'],
        default='gencode',
        help="Reference FASTA source (default: 'gencode')."
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        default='yes',
        required=False,
        help="If 'yes', gzip the output TSV file (default: 'yes')."
    )
    return parser.parse_args()


def run():
    # Step 1. Read the input arguments
    args = parse_args()

    # Step 2. Read the prob file
    transcripts = {} # key = index, value = transcript ID
    data = {
        'read_id': [],
        'transcript_id': [],
        'probability': []
    }
    with open(args.prob_file, "r") as f:
        num_transcripts = 0
        num_transcripts_iterated = 0
        first_line = True
        for line in f.readlines():
            if first_line:
                num_transcripts = int(line.split('\t')[0])
                first_line = False
            else:
                if num_transcripts_iterated < num_transcripts:
                    # Transcript ID
                    transcript_index = num_transcripts_iterated
                    if args.source == 'gencode':
                        transcript_id = line.replace('\n','')
                        transcript_id = transcript_id.split('|')[0]
                    else:
                        raise Exception('Unsupported source: %s', args.source)
                    transcripts[transcript_index] = transcript_id
                    num_transcripts_iterated += 1
                else:
                    # Read ID
                    elements = line.split('\t')
                    curr_read_id = str(elements[0])
                    curr_num_transcripts = int(elements[1])
                    curr_transcript_ids = []
                    curr_probabilities = []
                    for i in range(0, curr_num_transcripts):
                        curr_transcript_ids.append(int(elements[2+i]))
                        curr_probabilities.append(float(elements[2+curr_num_transcripts+i]))
                    for i,curr_transcript_id in enumerate(curr_transcript_ids):
                        data['read_id'].append(curr_read_id)
                        data['transcript_id'].append(transcripts[curr_transcript_ids[i]])
                        data['probability'].append(curr_probabilities[i])

    df = pd.DataFrame(data)

    # Step 3. Apply probability threshold
    df = df[df['probability'] >= args.threshold]

    # Step 4. Write the DataFrame to a TSV file
    if args.gzip:
        if args.output_tsv_file.endswith('.gz'):
            output_tsv_file = args.output_tsv_file
        else:
            output_tsv_file = args.output_tsv_file + '.gz'
    else:
        output_tsv_file = args.output_tsv_file

    df.to_csv(
        output_tsv_file,
        sep='\t',
        index=False,
        compression='gzip' if args.gzip else None
    )
