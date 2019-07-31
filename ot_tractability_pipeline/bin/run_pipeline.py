from ot_tractability_pipeline.buckets import *
import pandas as pd
import argparse
import datetime

def run(ensembl_id_list, database_url, out_file_name):

    #smab_buckets = pd.read_csv('tractability_buckets.tsv', sep='\t')
    # Assign tractability buckets

    # From Ensembl gene ID, get Uniprot, GO.CC, Symbol and PDB codes
    setup = Pipeline_setup(ensembl_id_list)

    # Get small molecule tractability info
    sm = Small_molecule_buckets(setup, database_url=database_url)
    sm_out_buckets = sm.assign_buckets()
    print(sm_out_buckets.groupby('Top_bucket')['ensembl_gene_id'].count())

    # Get antibody tractability info
    ab = Antibody_buckets(setup,database_url=database_url,sm_output=sm_out_buckets)
    smab_buckets = ab.assign_buckets()
    print(smab_buckets.groupby('Top_bucket_ab')['accession'].count())

    # Get protac tractability info
    protac = Protac_buckets(setup, ab_output=smab_buckets)
    out_buckets = protac.assign_buckets()
    print(out_buckets.groupby('Top_bucket_PROTAC')['accession'].count())


    d= datetime.date.today()
    out_buckets.to_csv('tractability_buckets_{}.tsv'.format(d), sep='\t')


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Assess tractability of gene list')
    parser.add_argument('genes',
                        help='A file of Ensembl Gene IDs, one per line and no header')
    parser.add_argument('--db',
                        help='Address to your local ChEMBL installation. Preferred method is to set the CHEMBL_DB environment variable',
                        default=None)

    parser.add_argument('--out_file', default='tractability_buckets.tsv',
                        help='Name of output csv file')

    args = parser.parse_args()

    # Get a unique list of Gene IDs from OT
    gene_list = pd.read_csv(args.genes, encoding='utf-8', header=None, names=['ensembl_gene_id'], usecols=[0])
    ensembl_id_list = list(gene_list['ensembl_gene_id'].unique())

    # URL to local ChEMBL database
    database_url = args.db

    run(ensembl_id_list, database_url, out_file_name=args.out_file)


if __name__ == '__main__':
    main()
