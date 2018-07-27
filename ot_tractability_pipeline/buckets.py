import json
import sys
import time
import zipfile
import io
import pkg_resources
import os
import re

import mygene
import pandas as pd
from sqlalchemy import create_engine

from ot_tractability_pipeline.sm_queries import *
from ot_tractability_pipeline.ab_queries import *

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2

DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline', 'data/')


class Pipeline_setup(object):
    '''
    Information retrieved from MyGene is used in both the small molecule and antibody pipelines. This class handles
    the aggregation of data ahead of assigning targets to buckets
    '''

    def __init__(self, ensembl_gene_id_list):

        # list of ensembl IDs for targets to be considered
        self.gene_list = ensembl_gene_id_list
        self._add_uniprot_column()

    def _uniprot_primary_only(self, a):
        '''
        If multiple uniprot IDs, only take primary (assume first in list is primary) <== NEEDS CHECKING
        :return:
        '''
        if isinstance(a, dict):
            s = a['Swiss-Prot']
            if isinstance(s, list):
                return s[0]
            else:
                return s
        else:
            return a

    def _add_uniprot_column(self):
        '''
        Use MyGene to find uniprot accession, entrezgene  pdb, pfam, GO and interpro IDs associated to each Ensembl ID.
        :return:
        '''

        mg = mygene.MyGeneInfo()

        # Use MyGene to return list of Uniprot accession numbers for each ensemble gene ID

        results = mg.getgenes(list(self.gene_list), scopes='ensembl',
                              as_dataframe=True, fields='uniprot.Swiss-Prot,entrezgene,pdb,go',
                              species='human', returnall=True)

        self.id_xref = results
        self.id_xref['uniprot'] = self.id_xref['uniprot'].apply(self._uniprot_primary_only)
        self.id_xref.reset_index(inplace=True)
        # self.id_xref.rename({'_id', '_score', 'entrezgene', 'go', 'interpro', 'pdb', 'pfam','uniprot'})
        self.id_xref.rename(columns={'query': 'ensembl_gene_id', 'uniprot': 'accession'}, inplace=True)


class Small_molecule_buckets(object):
    '''
    Class for assigning genes to tractability buckets
    '''

    ##############################################################################################################
    #
    # Initial setup
    #
    #
    ##############################################################################################################

    def __init__(self, Pipeline_setup, database_url=None, ligand_filter=list()):

        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # Cross referencing from Pipeline_setup, prevents repetition for antibody ot_tractability_pipeline

        self.id_xref = Pipeline_setup.id_xref

        # Load lists for PDB ligand filters (i.e. to remove sugars, solvents etc)

        self.ligand_filter = ligand_filter

        if len(ligand_filter) == 0:
            # Default ligand filter
            # Load unwanted PDB ligands list from text files

            with open(os.path.join(DATA_PATH, 'organic_solvents.txt')) as solvent_file:
                solvents = [a.split('\t')[0] for a in solvent_file]

            with open(os.path.join(DATA_PATH, 'sugars.txt')) as sugar_file:
                sugars = [a.split('\t')[0] for a in sugar_file]

            with open(os.path.join(DATA_PATH, 'cofactors.txt')) as cofactor_file:
                cofactors = [a.split('\t')[0] for a in cofactor_file]

            ligand_filter = solvents + sugars + cofactors

            self.ligand_filter = ligand_filter

        # All chembl data loaded into here
        self.all_chembl_targets = None

        # Unique list of PDB codes:
        self.pdb_list = []

        # Map back to Uniprot accession
        self.pdb_map = {}

        # If database url not supplied, get from envronemnt variable

        if database_url is None:
            database_url = os.getenv('CHEMBL_DB')
            print(database_url)

        # Create ChEMBL DB connection

        try:
            self.engine = create_engine(database_url)
        except AttributeError:
            print('''
Please supply a valid database URL to your local ChEMBL installation using one of the following methods

1) Supply the URL to this class within your script
2) Create an environment variable named 'CHEMBL_DB' with containing your database URL
3) Use the command line flag '--db' if using the supplied run script  '''
                  )
            raise

        # URL for PDBe web services

        self.PDB_SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    #
    #
    ##############################################################################################################

    def _search_chembl_clinical(self):
        '''
        Search for all targets in ChEMBL, returning data required for buckets 1-3
        :return:
        '''

        small_mol_info = pd.read_sql_query(chembl_clinical_small_mol, self.engine)
        self.all_chembl_targets = pd.read_sql_query(chembl_clinical_targets, self.engine)
        self.all_chembl_targets = self.all_chembl_targets.merge(small_mol_info, on='parent_molregno')

    def _process_protein_complexes(self):
        '''
        For protein complexes, see if we know the binding subunit, and only keep these

        :return:
        '''

        pc = self.all_chembl_targets[self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]
        not_pc = self.all_chembl_targets[~self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]

        n = 1000
        targets = pc['tid'].unique()
        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]

        df_list = []

        # Check if binding sites are defined
        for chunk in chunks:
            q = '''
            select distinct bs.site_id, td.tid 
            from CHEMBL_24.target_dictionary td, CHEMBL_24.binding_sites bs
            where td.tid = bs.tid and td.tid IN {}
            
            '''.format(tuple(chunk))
            df_list.append(pd.read_sql_query(q, self.engine))

        # Merge will set those with unknown binding site as NAN
        binding_site_info = pd.concat(df_list, sort=False)
        pc = pc.merge(binding_site_info, how='left', on='tid')
        defined = pc[pc['site_id'].notnull()]
        undefined = pc[~pc['site_id'].notnull()]

        # if binding site is defined, only take the subunits that are involved in the binding

        n = 1000
        targets = defined['accession'].unique()
        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]
        df_list2 = []

        for chunk in chunks:
            q2 = '''
            select distinct sc.component_id, cs.accession
            from CHEMBL_24.component_sequences cs, CHEMBL_24.site_components sc
            where cs.component_id = sc.component_id
            and cs.accession in {}'''.format(tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)
        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc], sort=False)

    def _assign_buckets_1_to_3(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        Possible duplication of effort of the OT known_drug score

        :return:
        '''

        def other_func(x):
            return tuple(set(x))

        self._search_chembl_clinical()
        self._process_protein_complexes()

        self.gene_xref = self.id_xref[['accession', 'ensembl_gene_id']]

        self.out_df = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')

        self.out_df.drop(['component_id', 'drug_name', 'ref_id', 'ref_type', 'tid', 'molregno',
                          'parent_molregno', 'ref_url'], axis=1, inplace=True)

        f = {x: 'first' for x in self.out_df.columns}
        f['max_phase'] = 'max'
        f['pref_name'] = other_func
        f['moa_chembl'] = other_func

        self.out_df = self.out_df.groupby(['ensembl_gene_id']).agg(f).reset_index(drop=True)

        self.out_df['Bucket_1'] = 0
        self.out_df['Bucket_2'] = 0
        self.out_df['Bucket_3'] = 0

        self.out_df.loc[(self.out_df['max_phase'] == 4), 'Bucket_1'] = 1
        self.out_df.loc[(self.out_df['max_phase'] < 4) & (self.out_df['max_phase'] >= 2), 'Bucket_2'] = 1
        self.out_df.loc[(self.out_df['max_phase'] < 2) & (self.out_df['max_phase'] > 0), 'Bucket_3'] = 1

    ##############################################################################################################
    #
    # Functions relating to bucket 4
    #
    #
    ##############################################################################################################

    def _make_request(self, url, data):
        request = urllib2.Request(url)

        try:
            url_file = urllib2.urlopen(request, data)
        except urllib2.HTTPError as e:
            if e.code == 404:
                print("[NOTFOUND %d] %s" % (e.code, url))
            else:
                print("[ERROR %d] %s" % (e.code, url))

            return None

        return url_file.read().decode()

    def _post_request(self, url, data, pretty=False):
        full_url = "%s/%s/?pretty=%s" % (self.PDB_SERVER_URL, url, str(pretty).lower())

        if isinstance(data, (list, tuple)):
            data = ",".join(data)

        return self._make_request(full_url, data.encode())

    def _pdb_list(self, s):

        pdb = s['pdb']
        acc = s['accession']

        if isinstance(pdb, list):
            self.pdb_list += pdb
            for p in pdb:
                p = p.lower()
                self.pdb_map[p] = acc
        elif isinstance(pdb, str):
            pdb = pdb.lower()
            self.pdb_list.append(pdb)
            self.pdb_map[pdb] = acc

    def _has_ligands(self, ligand_li):

        ligands = [l['chem_comp_id'] for l in ligand_li if l['chem_comp_id'] not in self.ligand_filter
                   and 'ION' not in l['chem_comp_name']]

        if len(ligands) > 0:
            return True
        else:
            return False

    def _pdb_ligand_info(self):
        '''
        Use PDBe webservices to get ligand information. POST requests allow up to 1000 PDBs to be submitted each request
        :return:
        '''

        self.id_xref.apply(self._pdb_list, axis=1)
        unique_list = [p.lower() for p in set(self.pdb_list)]
        n = 1000
        chunks = [unique_list[i:i + n] for i in range(0, len(unique_list), n)]

        self.no_ligands = []
        self.good_ligands = []
        self.bad_ligands = []

        for chunk in chunks:
            ligand_url = '/pdb/entry/ligand_monomers'
            data = ','.join(chunk)

            results = json.loads(self._post_request(ligand_url, data, False))

            # PDBs without ligands are not returned
            chunk_no_ligand = [p for p in chunk if p not in results.keys()]

            # Separate PDBs that do or don't contain suitable ligands
            chunk_good_ligand = [p for p in results.keys() if self._has_ligands(results[p])]
            chunk_bad_ligand = [p for p in results.keys() if not self._has_ligands(results[p])]

            # Add chunk info

            self.no_ligands += chunk_no_ligand
            self.good_ligands += chunk_good_ligand
            self.bad_ligands += chunk_bad_ligand

            time.sleep(1.5)

    def _known_pdb_ligand(self, s):

        if s in self.acc_known_lig:
            return 1
        else:
            return 0

    def _assign_bucket_4(self):
        '''
        Does the target have a ligand-bound protein crystal structure?
        '''

        # Download ligand info from pdb
        self._pdb_ligand_info()

        # Accession numbers with PDB ligand
        self.acc_known_lig = list(set([self.pdb_map[p] for p in self.good_ligands]))

        self.out_df['Bucket_4'] = self.out_df['accession'].apply(self._known_pdb_ligand)

    ##############################################################################################################
    #
    # Functions relating to buckets 5-6
    #
    #
    ##############################################################################################################

    def _assign_bucket_5_and_6(self):
        '''
        Does the target have a DrugEBIlity ensemble score >=0.7 (bucket 5) or  0<score<0.7 (bucket 6)
        '''
        df = pd.read_csv(os.path.join(DATA_PATH, 'drugebility_scores.csv'))

        df = df.merge(self.gene_xref, on='accession', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max()
        df['ensemble'].fillna(-1, inplace=True)

        print(self.out_df)
        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id')
        self.out_df['Bucket_5'] = 0
        self.out_df['Bucket_6'] = 0

        self.out_df.loc[(self.out_df['ensemble'] >= 0.7), 'Bucket_5'] = 1
        self.out_df.loc[(self.out_df['ensemble'] > 0) & (self.out_df['ensemble'] < 0.7), 'Bucket_6'] = 1

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    #
    #
    ##############################################################################################################

    def _search_chembl(self):

        self.activities = pd.concat([pd.read_sql_query(pchembl_q, self.engine),
                                     pd.read_sql_query(nm_q, self.engine),
                                     pd.read_sql_query(km_kon_q, self.engine),
                                     pd.read_sql_query(D_Tm_q, self.engine),
                                     pd.read_sql_query(residual_act_q, self.engine),
                                     pd.read_sql_query(Imax_q, self.engine),
                                     pd.read_sql_query(Emax_q, self.engine)], sort=False)

        # For faster testing
        # self.activities = pd.read_sql_query(pchembl_q, self.engine)

    def _structural_alerts(self):
        q = '''
        select distinct cs.canonical_smiles, csa.alert_id /*, sa.alert_name, sas.set_name */
        from CHEMBL_24.compound_structures cs,
            CHEMBL_24.compound_structural_alerts csa,
            CHEMBL_24.structural_alerts sa,
            CHEMBL_24.molecule_dictionary md,
            CHEMBL_24.structural_alert_sets sas
        where cs.molregno = md.molregno
        and cs.molregno = csa.molregno
        and csa.alert_id = sa.alert_id
        and sa.alert_set_id = 1
        and sa.alert_set_id = sas.alert_set_id
        '''

        alerts = pd.read_sql_query(q, self.engine)
        alerts = alerts.groupby('canonical_smiles', as_index=False).count()

        return alerts

    def _calc_pfi(self, s):

        ar = s['aromatic_rings']
        logd = s['acd_logd']
        return ar + logd

    def _assign_bucket_7(self):
        '''
        Does the target have ligands in ChEMBL (PFI <=7, SMART hits <= 2, scaffolds >= 2)
        '''
        self._search_chembl()
        self.activities['pfi'] = self.activities.apply(self._calc_pfi, axis=1)

        alerts = self._structural_alerts()

        self.activities = self.activities.merge(alerts, how='left', on='canonical_smiles')
        self.activities['alert_id'].fillna(0, inplace=True)
        self.activities.to_csv('activities.csv')
        df = self.activities[(self.activities['pfi'] <= 7) & (self.activities['alert_id'] <= 2)]

        f = {x: 'first' for x in df.columns}
        f['canonical_smiles'] = 'count'

        df2 = df.groupby('accession').agg(f).reset_index(drop=True)
        df2 = df2[['accession', 'canonical_smiles', 'target_chembl_id']]
        self.out_df = df2.merge(self.out_df, how='right', on='accession')
        self.out_df['Bucket_7'] = 0
        self.out_df['target_chembl_id_y'] = self.out_df['target_chembl_id_y'].fillna(self.out_df['target_chembl_id_x'])
        self.out_df.rename(columns={'target_chembl_query_y': 'target_chembl_query'}, inplace=True)

        self.out_df.loc[(self.out_df['canonical_smiles'] >= 2), 'Bucket_7'] = 1

        # PandasTools.AddMoleculeColumnToFrame(self.activities,'canonical_smiles','molecule')
        # PandasTools.AddMurckoToFrame(self.activities,molCol='molecule',MurckoCol='scaffold',Generic=True)
        # self.activities.to_csv('scaffolds.csv')

    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Is this target considered druggable using Finan et al's druggable genome?
    #
    ##############################################################################################################

    def _assign_bucket_8(self):
        '''
        Is this target considered druggable using Finan et al's druggable genome?
        '''
        df = pd.read_csv(os.path.join(DATA_PATH, 'druggable_genome.csv'))
        df = df[['ensembl_gene_id', 'small_mol_druggable']]
        df['small_mol_druggable'].fillna('N', inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id')
        self.out_df['Bucket_8'] = 0
        self.out_df.loc[(self.out_df['small_mol_druggable'] == 'Y'), 'Bucket_8'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 9
    #
    #
    ##############################################################################################################

    def _assign_bucket_9(self):
        '''
        Search to see whether targets have chemical patents in last 5 years

        '''

        pass

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _summarise_buckets(self):

        self.out_df['Bucket_sum'] = self.out_df['Bucket_1'] + self.out_df['Bucket_2'] + self.out_df[
            'Bucket_3'] + self.out_df['Bucket_4'] + self.out_df['Bucket_5'] + self.out_df['Bucket_6'] + self.out_df[
                                        'Bucket_7'] + self.out_df['Bucket_8']

        self.out_df['Top_bucket'] = 10
        for x in range(8, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}'.format(x)] == 1), 'Top_bucket'] = x

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._assign_buckets_1_to_3()
        self._assign_bucket_4()
        self._assign_bucket_5_and_6()
        self._assign_bucket_7()
        self._assign_bucket_8()
        # self._assign_bucket_9()
        self._summarise_buckets()

        self.out_df = self.out_df[['ensembl_gene_id', 'accession',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4', 'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_1', 'Bucket_sum', 'Top_bucket',
                                   'ensemble', 'canonical_smiles', 'small_mol_druggable']]

        self.out_df['Category'] = 'Unknown'

        self.out_df.loc[(self.out_df['Top_bucket'] <= 3), 'Category'] = 'Clinical Precedence'
        self.out_df.loc[(self.out_df['Top_bucket'] == 4) | (self.out_df['Top_bucket'] == 7),
                        'Category'] = 'Discovery Precedence'

        self.out_df.loc[
            (self.out_df['Top_bucket'] == 5) | (self.out_df['Top_bucket'] == 6) | (self.out_df['Top_bucket'] == 8),
            'Category'] = 'Predicted Tractable'

        return self.out_df


class Antibody_buckets(object):
    '''
    Class for assigning genes to tractability buckets
    '''

    ##############################################################################################################
    #
    # Initial setup
    #
    #
    ##############################################################################################################

    def __init__(self, Pipeline_setup, database_url=None, sm_output=None):

        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # Cross referencing from Pipeline_setup, prevents repetition for antibody ot_tractability_pipeline

        self.id_xref = Pipeline_setup.id_xref

        # If antibody results are to be combined with small molecule results, append antibody columns to sm results
        # Otherwise, use the id_xref dataframe

        # Load accepted GO locations
        self.accepted_go_locs = {}
        with open(os.path.join(DATA_PATH, 'go_accepted_loc.tsv')) as go_loc_file:
            for line in go_loc_file:
                line = line.split('\t')
                self.accepted_go_locs[line[0]] = line[2]

        if sm_output is not None:
            go_data = self.id_xref[['ensembl_gene_id', 'go']]
            self.out_df = sm_output.merge(go_data, how='outer', on='ensembl_gene_id')

        else:
            self.out_df = self.id_xref

        # All chembl data loaded into here
        self.all_chembl_targets = None

        # Unique list of PDB codes:
        self.pdb_list = []

        # Map back to Uniprot accession
        self.pdb_map = {}

        if database_url is None:
            database_url = os.getenv('CHEMBL_DB')
            print(database_url)

        # Create ChEMBL DB connection
        self.engine = create_engine(database_url)

    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    #
    #
    ##############################################################################################################

    def _process_protein_complexes(self):
        '''
        For protein complexes, see if we know the binding subunit, and only keep these

        :return:
        '''

        pc = self.all_chembl_targets[self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]
        not_pc = self.all_chembl_targets[~self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]

        defined = pc[pc['site_id'].notnull()]
        undefined = pc[~pc['site_id'].notnull()]

        # if binding site is defined, only take the subunits that are involved in the binding

        n = 1000
        targets = defined['site_id'].unique()
        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]
        df_list2 = []

        for chunk in chunks:
            q2 = '''
            select distinct cs.accession, sc.component_id
            from CHEMBL_24.component_sequences cs,
                CHEMBL_24.site_components sc
            where cs.component_id = sc.component_id
            and sc.site_id in {}'''.format(tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)
        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc], sort=False)

    def _assign_buckets_1_to_3(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        Possible duplication of effort of the OT known_drug score

        :return:
        '''

        self.all_chembl_targets = pd.read_sql_query(chembl_clinical_ab_targets, self.engine)
        self.all_chembl_targets.loc[self.all_chembl_targets['ref_type'] == 'Expert', ['ref_id', 'ref_url']] = 'NA'

        #

        self._process_protein_complexes()

        ab_info = pd.read_sql_query(chembl_clinical_ab, self.engine)
        self.all_chembl_targets = self.all_chembl_targets.merge(ab_info, how='left', on='parent_molregno')

        # Make sure max phase is for correct indication

        self.all_chembl_targets = self.all_chembl_targets[
            self.all_chembl_targets['max_phase'] == self.all_chembl_targets['max_phase_for_ind']]

        def other_func(x):
            return tuple(x)

        f = {x: other_func for x in self.all_chembl_targets if x != 'accession'}
        f['max_phase_for_ind'] = 'max'
        f['max_phase'] = 'max'

        self.all_chembl_targets = self.all_chembl_targets.groupby('accession', as_index=False).agg(f)

        self.out_df = self.all_chembl_targets.merge(self.out_df, how='outer', on='accession', suffixes=('_sm', '_ab'))

        self.out_df.drop(['component_id', 'drug_name', 'ref_id', 'ref_type', 'tid',
                          'ref_url'], axis=1, inplace=True)

        f = {x: 'first' for x in self.out_df.columns if x != 'ensembl_gene_id'}
        f['max_phase_for_ind'] = 'max'

        self.out_df = self.out_df.groupby(['ensembl_gene_id'], as_index=False).agg(f)

        self.out_df['Bucket_1_ab'] = 0
        self.out_df['Bucket_2_ab'] = 0
        self.out_df['Bucket_3_ab'] = 0

        self.out_df.loc[(self.out_df['max_phase_for_ind'] == 4), 'Bucket_1_ab'] = 1
        self.out_df.loc[
            (self.out_df['max_phase_for_ind'] < 4) & (self.out_df['max_phase_for_ind'] >= 2), 'Bucket_2_ab'] = 1
        self.out_df.loc[
            (self.out_df['max_phase_for_ind'] < 2) & (self.out_df['max_phase_for_ind'] > 0), 'Bucket_3_ab'] = 1

        self.out_df.to_csv('test.csv')

    ##############################################################################################################
    #
    # Functions relating to buckets 4, 6 and 7 (Uniprot loc)
    #
    #
    ##############################################################################################################

    def _make_request(self, url, data):
        request = urllib2.Request(url)

        try:
            url_file = urllib2.urlopen(request)
        except urllib2.HTTPError as e:
            if e.code == 404:
                print("[NOTFOUND %d] %s" % (e.code, url))
            else:
                print("[ERROR %d] %s" % (e.code, url))

            return None

        return url_file.read().decode()

    def _post_request(self, url, data):
        base = 'http://www.uniprot.org'
        full_url = "%s/%s" % (base, url)

        if isinstance(data, (list, tuple)):
            data = ",".join(data)

        return self._make_request(full_url, data.encode())

    def split_loc(self, s):

        try:
            loc_list = [a.split('. ') for a in s.split(';')]
            loc_evidence_li = []
            for l in loc_list:
                for l2 in l:
                    if l2 == '':
                        continue

                    l2 = l2.strip()
                    if l2.startswith('Note'):
                        break

                    l2 = l2.replace('SUBCELLULAR LOCATION: ', '')
                    #evidence = l2[l2.find("{") + 1:l2.find("}")]
                    evidence = re.findall(r'\{([^]]*)\}', l2)
                    print(evidence)
                    if len(evidence) == 0:
                        evidence = ['Unknown evidence type']


                    locations = l2.split('{')[0]

                    if not locations.startswith('Note') and not locations =='':
                        loc_evidence_li.append((evidence, locations))
            return loc_evidence_li

        except AttributeError:
            return [('na', 'na')]

    def _set_b4_flag(self, s):


        accepted_uniprot_high_conf = [a[1] for a in s['Subcellular location [CC]'] if
                                      ('Cell membrane' in a[1] or 'Secreted' in a[1]) and (
                                                  'ECO:0000269' in a[0] or 'ECO:0000305' in a[0])]

        all_uniprot_high_conf = {a[1]:a[0] for a in s['Subcellular location [CC]'] if
                                 ('ECO:0000269' in a[0] or 'ECO:0000305' in a[0])}

        if len(accepted_uniprot_high_conf) > 0:
            b4_flag = 1
        else:
            b4_flag = 0

        return b4_flag, all_uniprot_high_conf

    def _set_b6_flag(self, s):

        accepted_uniprot_med_conf = [a[1] for a in s['Subcellular location [CC]'] if
                                     ('Cell membrane' in a[1] or 'Secreted' in a[1]) and not (
                                             'ECO:0000269' in a[0] or 'ECO:0000305' in a[0])]

        all_uniprot_med_conf = {a[1]:a[0] for a in s['Subcellular location [CC]'] if not
        ('ECO:0000269' in a[0] or 'ECO:0000305' in a[0])}

        if len(accepted_uniprot_med_conf) > 0:
            b6_flag = 1
        else:
            b6_flag = 0

        return b6_flag, all_uniprot_med_conf

    def _assign_bucket_4_and_6(self):
        '''
        Uniprot (loc): Targets in "Cell membrane" or "Secreted", high confidence
        '''

        url = "uniprot/?format=tab&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&columns=id,comment(SUBCELLULAR+LOCATION),comment(DOMAIN),feature(DOMAIN+EXTENT),feature(INTRAMEMBRANE),feature(TOPOLOGICAL+DOMAIN),feature(TRANSMEMBRANE),feature(SIGNAL)"
        data = ['P42336', 'P60484']

        location = self._post_request(url, data)
        location = [x.split('\t') for x in location.split('\n')]
        df = pd.DataFrame(location[1:], columns=location[0])
        df['uniprot_loc_test']= df['Subcellular location [CC]']

        df['Subcellular location [CC]'] = df['Subcellular location [CC]'].apply(self.split_loc)
        df['Bucket_4_ab'], df['Uniprot_high_conf_loc'] = zip(*df.apply(self._set_b4_flag, axis=1))
        df['Bucket_6_ab'], df['Uniprot_med_conf_loc'] = zip(*df.apply(self._set_b6_flag, axis=1))
        df.rename(columns={'Entry': 'accession'}, inplace=True)

        self.out_df = self.out_df.merge(df, how='left', on='accession')

    ##############################################################################################################
    #
    # Functions relating to buckets 5 and 8 (GO CC)
    #
    #
    ##############################################################################################################

    def _set_b5_b8_flag(self, s):
        try:
            cc = s['go']['CC']
        except:
            return 0, [], 0, []

        evidence_types = {'EXP': 'High', 'IDA': 'High', 'IPI': 'High', 'TAS': 'High', 'IMP': 'High', 'IGI': 'High',
                          'IEP': 'High',
                          'ISS': 'Medium', 'ISO': 'Medium', 'ISA': 'Medium', 'ISM': 'Medium', 'IGC': 'Medium',
                          'IBA': 'Medium', 'IBD': 'Medium', 'IKR': 'Medium', 'IRD': 'Medium', 'RCA': 'Medium',
                          'IEA': 'Medium',
                          'NAS': 'Low', 'IC': 'Low', 'ND': 'Low', 'NR': 'Low'
                          }

        high_conf_loc = []
        med_conf_loc = []
        accepted_high_conf_loc = []
        accepted_med_conf_loc = []

        if isinstance(cc, dict):
            cc = [cc]

        for c_dict in cc:
            try:
                go_id = c_dict['id']
                go_loc = c_dict['term']
                evidence = c_dict['evidence']
            except TypeError:
                continue
            try:
                confidence = evidence_types[evidence]
            except KeyError:
                confidence = None

            if confidence == 'High':
                high_conf_loc.append((go_loc, evidence))
            elif confidence == 'Medium':
                med_conf_loc.append((go_loc, evidence))

            if go_id in self.accepted_go_locs.keys():
                if confidence == 'High':
                    accepted_high_conf_loc.append(go_loc)
                elif confidence == 'Medium':
                    accepted_med_conf_loc.append(go_loc)

        b5_flag = 0
        b8_flag = 0

        if len(accepted_high_conf_loc) > 0:
            b5_flag = 1
        elif len(accepted_med_conf_loc) > 0:
            b8_flag = 1

        return b5_flag, high_conf_loc, b8_flag, med_conf_loc

    def _assign_bucket_5_and_8(self):
        '''
        GO CC
        '''
        self.out_df['Bucket_5_ab'], self.out_df['GO_high_conf_loc'], self.out_df['Bucket_8_ab'], self.out_df[
            'GO_med_conf_loc'] = zip(*self.out_df.apply(self._set_b5_b8_flag, axis=1))

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    #
    #
    ##############################################################################################################

    def _assign_bucket_7(self):
        '''
        Uniprot (SigP + TMHMM): targets with predicted Signal Peptide or Trans-membrane regions, and not destined to
        organelles

        '''
        self.out_df['Bucket_7_ab'] = 0

        self.out_df.loc[(self.out_df['Transmembrane'].str.contains('TRANSMEM', na=False)), 'Bucket_7_ab'] = 1
        self.out_df.loc[(self.out_df['Signal peptide'].str.contains('SIGNAL', na=False)), 'Bucket_7_ab'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 9
    #
    #
    ##############################################################################################################

    def _main_location(self, s):

        if s['Reliability'] == 'Validated':
            return s['Validated']
        else:
            return s['Supported']

    def _assign_bucket_9(self):
        '''
        HPA
        '''

        with urllib2.urlopen(
                'https://www.proteinatlas.org/download/subcellular_location.tsv.zip') as zip_file:
            with zipfile.ZipFile(io.BytesIO(zip_file.read()), 'r') as pa_file:
                with pa_file.open('subcellular_location.tsv') as subcell_loc:
                    df = pd.read_csv(subcell_loc, sep='\t', header=0)

        df.rename(columns={'Gene': 'ensembl_gene_id'}, inplace=True)

        df['main_location'] = df.apply(self._main_location, axis=1)
        reliable = df[(df['Reliability'] == 'Supported') | (df['Reliability'] == 'Validated')]

        reliable.drop(columns=["Approved", "Supported", "Uncertain",
                               "Cell cycle dependency", "GO id"], axis=1, inplace=True)

        self.out_df = self.out_df.merge(reliable, on='ensembl_gene_id', how='left')

        self.out_df['Bucket_9_ab'] = 0
        self.out_df.loc[(self.out_df['main_location'].str.contains("Plasma membrane", na=False)), 'Bucket_9_ab'] = 1

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _summarise_buckets(self):

        self.out_df.drop('go', inplace=True, axis=1)

        self.out_df['Bucket_sum_ab'] = self.out_df['Bucket_1_ab'] + self.out_df['Bucket_2_ab'] + self.out_df[
            'Bucket_3_ab'] + self.out_df['Bucket_4_ab'] + self.out_df['Bucket_5_ab'] + self.out_df['Bucket_6_ab'
                                       ] + self.out_df['Bucket_7_ab'] + self.out_df['Bucket_8_ab'] + self.out_df[
                                           'Bucket_9_ab']

        self.out_df['Top_bucket_ab'] = 10
        for x in range(9, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_ab'.format(x)] == 1), 'Top_bucket_ab'] = x

        self.out_df.set_index('ensembl_gene_id')

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._assign_buckets_1_to_3()
        self._assign_bucket_4_and_6()
        self._assign_bucket_5_and_8()
        self._assign_bucket_7()
        self._assign_bucket_9()

        self._summarise_buckets()

        # try:

        self.out_df.index = self.out_df['ensembl_gene_id']

        self.out_df = self.out_df[['accession',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4',
                                   'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_sum', 'Top_bucket','Category',
                                   'ensemble', 'canonical_smiles', 'small_mol_druggable',
                                   'Bucket_1_ab', 'Bucket_2_ab', 'Bucket_3_ab', 'Bucket_4_ab',
                                   'Bucket_5_ab', 'Bucket_6_ab', 'Bucket_7_ab',
                                   'Bucket_8_ab', 'Bucket_9_ab', 'Bucket_sum_ab', 'Top_bucket_ab',
                                   'uniprot_loc_test','Uniprot_high_conf_loc', 'GO_high_conf_loc', 'Uniprot_med_conf_loc',
                                   'GO_med_conf_loc', 'Transmembrane', 'Signal peptide', 'main_location'
                                   ]]
        self.out_df.rename(columns={'canonical_smiles': 'High Quality ChEMBL compounds',
                                    'small_mol_druggable': 'Small Molecule Druggable Genome Member',
                                    'main_location': 'HPA main location'}, inplace=True)
        self.out_df.sort_values(['Top_bucket','Bucket_sum'], ascending=[True, False], inplace=True)



        return self.out_df
