from opentargets import OpenTargetsClient
from pipeline.buckets import *


disease= 'type 2 diabetes mellitus'

ot = OpenTargetsClient()
a_for_disease = ot.get_associations_for_disease(disease).to_dataframe()
a_for_disease.rename(columns={'target.id': 'ensembl_gene_id'}, inplace=True)

ensembl_id_list = list(a_for_disease['ensembl_gene_id'].unique())
setup = Pipeline_setup(ensembl_id_list)

setup.id_xref.to_csv('diabetes_xref.csv')