from pipeline.buckets import *
import pandas as pd

# Get a unique list of Gene IDs from OT
ot_data = pd.read_csv('all_targets.csv', encoding='utf-8')
ensembl_id_list = list(ot_data['ensembl_gene_id'].unique())

# URL to local ChEMBL database
database_url = 'oracle://cradoux:cradoux@ora-vm-065.ebi.ac.uk:1531/chempro'

# Assign tractability buckets
setup = Pipeline_setup(ensembl_id_list)

#sm = Small_molecule_buckets(setup, database_url=database_url)

sm = pd.read_csv('out_buckets2.csv')
ab = Antibody_buckets(setup,database_url=database_url)

sm_out_buckets = sm.assign_buckets()
#print(sm_out_buckets.groupby('Bucket')['ensembl_gene_id'].count())

#ab_out_buckets = ab.assign_buckets()
#print(ab_out_buckets.groupby('Bucket')['ensembl_gene_id'].count())

sm_out_buckets.to_csv('ab_out_buckets.csv')