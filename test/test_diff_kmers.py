from exactolib.main import diff_kmers


def test_diff_kmers_1():
    query_sequences = {
        'query_1': 'MSIINFEKLFG'
    }
    reference_sequences = {
        'ref_1': 'MAAKSDGRLKMKKSSDVAFTPLQNSDNSGSVQG',
        'ref_2': 'MVILLENTALSALVIIQFCWHFFATAFSLPPEE',
        'ref_3': 'MFVHDFSTEDSSTTAGSGGSAGSGGSKMVVIAA'
    }
    df_unique_kmers = diff_kmers(
        query_sequences=query_sequences,
        reference_sequences=reference_sequences,
        min_k=8,
        max_k=11
    )
    assert(len(df_unique_kmers) == 10)

