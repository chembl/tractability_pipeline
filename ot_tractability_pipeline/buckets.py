import io
import json
import re
import sys
import time
import zipfile

import mygene
import numpy as np
import pandas as pd
from pandas.io.json import json_normalize
import pkg_resources
from sqlalchemy import create_engine

from ot_tractability_pipeline.ab_queries import *
from ot_tractability_pipeline.sm_queries import *

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

    def __init__(self, ensembl_gene_id_list, store_fetched):

        # list of ensembl IDs for targets to be considered
        self.store_fetched = store_fetched
        self.gene_list = ensembl_gene_id_list
        self._add_uniprot_column()
        

    def _uniprot_primary_only(self, s):
        '''
        If multiple uniprot IDs, only take primary (assume first in list is primary) <== NEEDS CHECKING
        :return:
        '''

        if isinstance(s, list):
            return s[0]
        else:
            return s

    def _add_uniprot_column(self):
        '''
        Use MyGene to find uniprot accession, entrezgene  pdb, pfam, GO and interpro IDs associated to each Ensembl ID.
        :return:
        '''

        mg = mygene.MyGeneInfo()

        # Use MyGene to return list of Uniprot accession numbers for each ensemble gene ID

        results = mg.getgenes(list(self.gene_list), scopes='ensembl',
                              as_dataframe=True, fields='symbol,uniprot.Swiss-Prot,entrezgene,pdb,go',
                              species='human', returnall=True)

        self.id_xref = results

        self.id_xref['uniprot'] = self.id_xref['uniprot.Swiss-Prot'].apply(self._uniprot_primary_only)
        self.id_xref.reset_index(inplace=True)
        # self.id_xref.rename({'_id', '_score', 'entrezgene', 'go', 'interpro', 'pdb', 'pfam','uniprot'})
        self.id_xref.rename(columns={'query': 'ensembl_gene_id', 'uniprot': 'accession'}, inplace=True)
        print(self.id_xref.columns)

        if self.store_fetched: 
            self.id_xref.to_csv("{}/id_xref.csv".format(self.store_fetched))


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
        self.store_fetched = Pipeline_setup.store_fetched
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
        self.acc_map = {}

        # If database url not supplied, get from envronemnt variable

        if database_url is None:
            database_url = os.getenv('CHEMBL_DB')

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
    # Clinical Compounds
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

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/sm_all_chembl_targets.csv".format(self.store_fetched))

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
            from {0}.target_dictionary td, {0}.binding_sites bs
            where td.tid = bs.tid and td.tid IN {1}
            
            '''.format(CHEMBL_VERSION, tuple(chunk))
            df_list.append(pd.read_sql_query(q, self.engine))

        # Merge will set those with unknown binding site as NAN
        binding_site_info = pd.concat(df_list, sort=False)

        if self.store_fetched: 
            binding_site_info.to_csv("{}/sm_chembl_binding_site_info.csv".format(self.store_fetched))

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
            from {0}.component_sequences cs, {0}.site_components sc
            where cs.component_id = sc.component_id
            and cs.accession in {1}'''.format(CHEMBL_VERSION, tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)

        if self.store_fetched: 
            binding_subunit.to_csv("{}/sm_chembl_binding_subunit.csv".format(self.store_fetched))

        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc], sort=False)

    def _assess_clinical(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        '''

        def other_func(x):
            return tuple(set(x))

        print(self.id_xref['symbol'])
        self._search_chembl_clinical()
        self._process_protein_complexes()

        self.gene_xref = self.id_xref[['accession', 'ensembl_gene_id', 'symbol']]

        self.out_df = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')

        self.clinical_evidence = self.out_df

        self.clinical_evidence.to_csv('clinical_evidence.csv', index=False)

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

        print(self.out_df['symbol'])

    ##############################################################################################################
    #
    # Functions relating to bucket 4
    # PDB
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

        if not isinstance(pdb, list): pdb = [pdb]

        # Python 2/3 compatability
        try: pdb = [p.lower() for p in pdb if isinstance(p,(str,unicode))] #Python 2
        except: pdb = [p.lower() for p in pdb if isinstance(p,str)] #Python 3

        self.pdb_list += pdb
        for p in pdb:
            try: self.pdb_map[p].add(acc)
            except KeyError: self.pdb_map[p] = {acc}
            
            try: self.acc_map[acc].append(p)
            except KeyError: self.acc_map[acc] = [p]

        self.pdb_list = list(set(self.pdb_list))

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

        # Fails with n=1000, runs with n=750
        n = 750
        chunks = [self.pdb_list[i:i + n] for i in range(0, len(self.pdb_list), n)]

        self.no_ligands = []
        self.good_ligands = []
        self.bad_ligands = []

        all_results = {}
        for chunk in chunks:
            ligand_url = '/pdb/entry/ligand_monomers'

            data = ','.join(chunk)

            results = json.loads(self._post_request(ligand_url, data, False))
            all_results.update(results)
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

        if self.store_fetched: 
            json.dump(all_results, open("{}/sm_pdb_ligand_info.json".format(self.store_fetched), 'wt'))


    def _known_pdb_ligand(self, s):

        if s in self.acc_known_lig:
            return list(set([p for p in self.acc_map[s] if p in self.good_ligands]))
        else:
            return np.nan

    def _assess_pdb(self):
        '''
        Does the target have a ligand-bound protein crystal structure?
        '''

        # Download ligand info from pdb
        self._pdb_ligand_info()

        # Accession numbers with PDB ligand
        self.acc_known_lig = list({c for pdb in self.good_ligands for c in self.pdb_map[pdb]})

        self.out_df['PDB_Known_Ligand'] = self.out_df['accession'].apply(self._known_pdb_ligand)

        self.out_df['Bucket_4'] = 0

        self.out_df.loc[(self.out_df['PDB_Known_Ligand'].notna()), 'Bucket_4'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 5-6
    # DrugEBIlity
    #
    ##############################################################################################################

    def _assess_pockets(self):
        '''
        Does the target have a DrugEBIlity ensemble score >=0.7 (bucket 5) or  0<score<0.7 (bucket 6)
        '''
        df = pd.read_csv(os.path.join(DATA_PATH, 'drugebility_scores.csv'))

        df = df.merge(self.gene_xref, on='accession', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max()
        df['ensemble'].fillna(-1, inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id', suffixes=['_drop', ''])
        self.out_df['Bucket_5'] = 0
        self.out_df['Bucket_6'] = 0

        self.out_df.loc[(self.out_df['ensemble'] >= 0.7), 'Bucket_5'] = 1
        self.out_df.loc[(self.out_df['ensemble'] > 0) & (self.out_df['ensemble'] < 0.7), 'Bucket_6'] = 1

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    # ChEBML
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

        if self.store_fetched: 
            self.activities.to_csv("{}/sm_chembl_activities.csv".format(self.store_fetched))

        # For faster testing
        # self.activities = pd.read_sql_query(pchembl_q, self.engine)

    def _structural_alerts(self):
        q = '''
        select distinct cs.canonical_smiles, csa.alert_id /*, sa.alert_name, sas.set_name */
        from {0}.compound_structures cs,
            {0}.compound_structural_alerts csa,
            {0}.structural_alerts sa,
            {0}.molecule_dictionary md,
            {0}.structural_alert_sets sas
        where cs.molregno = md.molregno
        and cs.molregno = csa.molregno
        and csa.alert_id = sa.alert_id
        and sa.alert_set_id = 1
        and sa.alert_set_id = sas.alert_set_id
        '''.format(CHEMBL_VERSION)

        alerts = pd.read_sql_query(q, self.engine)

        if self.store_fetched: 
            alerts.to_csv("{}/sm_chembl_alerts.csv".format(self.store_fetched))

        alerts = alerts.groupby('canonical_smiles', as_index=False).count()

        return alerts

    def _calc_pfi(self, s):

        ar = s['aromatic_rings']
        logd = s['acd_logd']
        return ar + logd

    def _assess_chembl(self):
        '''
        Does the target have ligands in ChEMBL (PFI <=7, SMART hits <= 2, scaffolds >= 2)
        Scaffold counting currently not implemented

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

        # Use RDKit to count scaffolds
        # PandasTools.AddMoleculeColumnToFrame(self.activities,'canonical_smiles','molecule')
        # PandasTools.AddMurckoToFrame(self.activities,molCol='molecule',MurckoCol='scaffold',Generic=True)
        # self.activities.to_csv('scaffolds.csv')

    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Is this target considered druggable using Finan et al's druggable genome?
    #
    ##############################################################################################################

    def _assess_druggable_genome(self):
        '''
        Is this target considered druggable using Finan et al's druggable genome?
        '''
        df = pd.read_csv(os.path.join(DATA_PATH, 'druggable_genome.csv'))
        df = df[['ensembl_gene_id', 'small_mol_druggable']]

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id')
        self.out_df['Bucket_8'] = 0
        self.out_df.loc[(self.out_df['small_mol_druggable'] == 'Y'), 'Bucket_8'] = 1
        self.out_df['small_mol_druggable'].fillna('N', inplace=True)

    ##############################################################################################################
    #
    # Functions relating to buckets 9
    #
    #
    ##############################################################################################################

    def _assign_bucket_9(self):
        '''
        Future bucket

        '''

        pass

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _summarise_buckets(self):

        '''
        Calculate the best highest bucket for each target
        :return:
        '''

        self.out_df['Top_bucket'] = 9
        for x in range(8, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}'.format(x)] == 1), 'Top_bucket'] = x
            self.out_df['Bucket_{}'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum'] = self.out_df['Bucket_1'] + self.out_df['Bucket_2'] + self.out_df[
            'Bucket_3'] + self.out_df['Bucket_4'] + self.out_df['Bucket_5'] + self.out_df['Bucket_6'] + self.out_df[
                                        'Bucket_7'] + self.out_df['Bucket_8']

    def _clinical_precedence(self, s):
        return 1 * s['Bucket_1'] + 0.7 * s['Bucket_2'] + 0.2 * s['Bucket_3']

    def _discovery_precedence(self, s):
        return 0.7 * s['Bucket_4'] + 0.3 * s['Bucket_7']

    def _predicted_tractable(self, s):
        return 0.7 * s['Bucket_5'] + 0.3 * s['Bucket_6'] + 0.3 * s['Bucket_8']

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._assess_clinical()
        self._assess_pdb()
        self._assess_pockets()
        self._assess_chembl()
        self._assess_druggable_genome()
        # self._assign_bucket_9()
        self._summarise_buckets()

        print(self.out_df.columns)
        # Add extra buckets to the list below
        self.out_df = self.out_df[['ensembl_gene_id', 'symbol', 'accession',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4', 'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_sum', 'Top_bucket',
                                   'ensemble', 'canonical_smiles', 'small_mol_druggable', 'PDB_Known_Ligand']]

        # Calculate category scores and assign highest category to each target

        self.out_df['Category'] = 'Unknown'
        self.out_df['Clinical_Precedence'] = self.out_df.apply(self._clinical_precedence, axis=1)
        self.out_df['Discovery_Precedence'] = self.out_df.apply(self._discovery_precedence, axis=1)
        self.out_df['Predicted_Tractable'] = self.out_df.apply(self._predicted_tractable, axis=1)

        self.out_df.loc[(self.out_df['Top_bucket'] <= 3), 'Category'] = 'Clinical_Precedence'
        self.out_df.loc[(self.out_df['Top_bucket'] == 4) | (self.out_df['Top_bucket'] == 7),
                        'Category'] = 'Discovery_Precedence'

        self.out_df.loc[
            (self.out_df['Top_bucket'] == 5) | (self.out_df['Top_bucket'] == 6) | (self.out_df['Top_bucket'] == 8),
            'Category'] = 'Predicted_Tractable'

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
        self.store_fetched = Pipeline_setup.store_fetched
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
            go_data = self.id_xref[['ensembl_gene_id', 'go.CC']]
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

        # Create ChEMBL DB connection
        self.engine = create_engine(database_url)

    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    # Clinical antibodies
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
            from {0}.component_sequences cs,
                {0}.site_components sc
            where cs.component_id = sc.component_id
            and sc.site_id in {1}'''.format(CHEMBL_VERSION, tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)

        if self.store_fetched: 
            binding_subunit.to_csv("{}/ab_chembl_binding_subunit.csv".format(self.store_fetched))

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

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/ab_all_chembl_targets.csv".format(self.store_fetched))

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

    ##############################################################################################################
    #
    # Functions relating to buckets 4, 6 and 7
    # Uniprot location
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
        '''
        Process the delimited string returned from uniprot webservice call
        :param s:
        :return:
        '''
        try:
            loc_list = [a.split('. ') for a in s.split(';')]
            loc_evidence_li = []

        except AttributeError:
            return [('na', 'na')]

        for l in loc_list:
            for l2 in l:
                if l2 == '':
                    continue

                l2 = l2.strip()
                if l2.startswith('Note'):
                    break

                l2 = l2.replace('SUBCELLULAR LOCATION: ', '')
                # evidence = l2[l2.find("{") + 1:l2.find("}")]
                evidence = re.findall(r'\{([^]]*)\}', l2)
                evidence = [e for x in evidence for e in x.split(',')]

                if len(evidence) == 0:
                    evidence = ['Unknown evidence type']

                locations = l2.split('{')[0]

                if locations != '':
                    loc_evidence_li.append((evidence, locations))
        return loc_evidence_li

    def _check_evidence(self, evidence_li):

        high_conf_evidence = [e for e in evidence_li if ('ECO:0000269' in e or 'ECO:0000305' in e)]

        if len(high_conf_evidence) > 0:
            return True
        return False

    def _set_b4_flag(self, s):

        accepted_uniprot_high_conf = [a[1] for a in s['Subcellular location [CC]'] if
                                      ('Cell membrane' in a[1] or 'Secreted' in a[1]) and (self._check_evidence(a[0]))]

        all_uniprot_high_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if self._check_evidence(a[0])]

        if len(accepted_uniprot_high_conf) > 0:
            b4_flag = 1
        else:
            b4_flag = 0

        return b4_flag, all_uniprot_high_conf

    def _set_b6_flag(self, s):

        accepted_uniprot_med_conf = [a[1] for a in s['Subcellular location [CC]'] if
                                     ('Cell membrane' in a[1] or 'Secreted' in a[1]) and not self._check_evidence(a[0])]

        all_uniprot_med_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if not self._check_evidence(a[0])]

        if len(accepted_uniprot_med_conf) > 0:
            b6_flag = 1
        else:
            b6_flag = 0

        return b6_flag, all_uniprot_med_conf

    def _assign_bucket_4_and_6(self):
        '''
        Uniprot (loc): Targets in "Cell membrane" or "Secreted", high confidence
        '''

        # Return all reviewed and Human targets
        url = "uniprot/?format=tab&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&columns=id,comment(SUBCELLULAR+LOCATION),comment(DOMAIN),feature(DOMAIN+EXTENT),feature(INTRAMEMBRANE),feature(TOPOLOGICAL+DOMAIN),feature(TRANSMEMBRANE),feature(SIGNAL)"
        data = ['P42336', 'P60484']

        location = self._post_request(url, data)
        location = [x.split('\t') for x in location.split('\n')]
        df = pd.DataFrame(location[1:], columns=location[0])
        df['uniprot_loc_test'] = df['Subcellular location [CC]']

        df['Subcellular location [CC]'] = df['Subcellular location [CC]'].apply(self.split_loc)

        if self.store_fetched: 
            df.to_csv("{}/ab_uniprot_for_buckets_4_and_6.csv".format(self.store_fetched))

        df['Bucket_4_ab'], df['Uniprot_high_conf_loc'] = zip(*df.apply(self._set_b4_flag, axis=1))
        df['Bucket_6_ab'], df['Uniprot_med_conf_loc'] = zip(*df.apply(self._set_b6_flag, axis=1))
        df.rename(columns={'Entry': 'accession'}, inplace=True)

        self.out_df = self.out_df.merge(df, how='left', on='accession')

    ##############################################################################################################
    #
    # Functions relating to buckets 5 and 8
    # GO Cell Component
    #
    ##############################################################################################################

    def _set_b5_b8_flag(self, s):
        try:
            cc = s['go.CC']
        except:
            return 0, [], 0, []

        # Confidence for each evidence type
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

        if not isinstance(cc, list):
            return 0, [], 0, []

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
    # Uniprot transmembrane and signal peptides
    #
    ##############################################################################################################

    def _split_loc_b7(self, s):

        try:
            loc_list = [a.split('.') for a in s.split(';')]
            loc_evidence_li = []

        except AttributeError:
            return [('na', 'na')]

        for l in loc_list:
            for l2 in l:
                if l2 == '':
                    continue

                l2 = l2.strip()
                l2 = re.sub(r'{[^}]+}', '', l2)
                if l2.startswith('TRANSMEM') or l2.startswith('SIGNAL'):
                    loc_evidence_li.append(l2)
        return loc_evidence_li

    def _assign_bucket_7(self):
        '''
        Uniprot (SigP + TMHMM): targets with predicted Signal Peptide or Trans-membrane regions, and not destined to
        organelles

        '''
        self.out_df['Bucket_7_ab'] = 0

        self.out_df.loc[(self.out_df['Transmembrane'].str.contains('TRANSMEM', na=False)), 'Bucket_7_ab'] = 1
        self.out_df.loc[(self.out_df['Signal peptide'].str.contains('SIGNAL', na=False)), 'Bucket_7_ab'] = 1

        self.out_df['Transmembrane'] = self.out_df['Transmembrane'].apply(self._split_loc_b7)
        self.out_df['Signal peptide'] = self.out_df['Signal peptide'].apply(self._split_loc_b7)

    ##############################################################################################################
    #
    # Functions relating to buckets 9
    # Human protein atlas - Main location
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

        # Download latest file
        zip_file = urllib2.urlopen('https://www.proteinatlas.org/download/subcellular_location.tsv.zip')
        with zipfile.ZipFile(io.BytesIO(zip_file.read()), 'r') as pa_file:
            with pa_file.open('subcellular_location.tsv') as subcell_loc:
                df = pd.read_csv(subcell_loc, sep='\t', header=0)

        df.rename(columns={'Gene': 'ensembl_gene_id'}, inplace=True)

        if self.store_fetched: 
            df.to_csv("{}/ab_proteinatlas_for_bucket_9.csv".format(self.store_fetched))

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

    def _clinical_precedence(self, s):
        return 1 * s['Bucket_1_ab'] + 0.7 * s['Bucket_2_ab'] + 0.2 * s['Bucket_3_ab']

    def _high_conf_pred(self, s):
        return 0.7 * s['Bucket_4_ab'] + 0.3 * s['Bucket_5_ab']

    def _med_conf_pred(self, s):
        return 0.4 * s['Bucket_6_ab'] + 0.25 * s['Bucket_7_ab'] + 0.25 * s['Bucket_8_ab'] + 0.1 * s['Bucket_9_ab']

    def _summarise_buckets(self):

        self.out_df.drop('go.CC', inplace=True, axis=1)

        self.out_df['Top_bucket_ab'] = 10
        for x in range(9, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_ab'.format(x)] == 1), 'Top_bucket_ab'] = x
            self.out_df['Bucket_{}_ab'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_ab'] = self.out_df['Bucket_1_ab'] + self.out_df['Bucket_2_ab'] + self.out_df[
            'Bucket_3_ab'] + self.out_df['Bucket_4_ab'] + self.out_df['Bucket_5_ab'] + self.out_df['Bucket_6_ab'
                                       ] + self.out_df['Bucket_7_ab'] + self.out_df['Bucket_8_ab'] + self.out_df[
                                           'Bucket_9_ab']

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

        # Columns to keep. This includes columns from the small molecule pipeline

        self.out_df = self.out_df[['symbol', 'accession',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4',
                                   'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_sum', 'Top_bucket', 'Category',
                                   'Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable',
                                   'PDB_Known_Ligand',
                                   'ensemble', 'canonical_smiles', 'small_mol_druggable',
                                   'Bucket_1_ab', 'Bucket_2_ab', 'Bucket_3_ab', 'Bucket_4_ab',
                                   'Bucket_5_ab', 'Bucket_6_ab', 'Bucket_7_ab',
                                   'Bucket_8_ab', 'Bucket_9_ab', 'Bucket_sum_ab', 'Top_bucket_ab',
                                   'Uniprot_high_conf_loc', 'GO_high_conf_loc',
                                   'Uniprot_med_conf_loc',
                                   'GO_med_conf_loc', 'Transmembrane', 'Signal peptide', 'main_location'
                                   ]]
        self.out_df.rename(columns={'canonical_smiles': 'High_Quality_ChEMBL_compounds',
                                    'small_mol_druggable': 'Small_Molecule_Druggable_Genome_Member',
                                    'main_location': 'HPA_main_location', 'Signal peptide': 'Signal_peptide'},
                           inplace=True)
        self.out_df.sort_values(['Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable'],
                                ascending=[False, False, False], inplace=True)

        # Score each category, and label highest category
        self.out_df['Clinical_Precedence_ab'] = self.out_df.apply(self._clinical_precedence, axis=1)
        self.out_df['Predicted_Tractable__High_confidence'] = self.out_df.apply(self._high_conf_pred, axis=1)
        self.out_df['Predicted_Tractable__Medium_to_low_confidence'] = self.out_df.apply(self._med_conf_pred, axis=1)

        self.out_df['Category_ab'] = 'Unknown'

        self.out_df.loc[(self.out_df['Top_bucket_ab'] <= 3), 'Category_ab'] = 'Clinical_Precedence'
        self.out_df.loc[(self.out_df['Top_bucket_ab'] == 4) | (self.out_df['Top_bucket_ab'] == 5),
                        'Category_ab'] = 'Predicted_Tractable__High_confidence'

        self.out_df.loc[
            (self.out_df['Top_bucket_ab'] == 6) | (self.out_df['Top_bucket_ab'] == 7) | (
                    self.out_df['Top_bucket_ab'] == 8) | (self.out_df['Top_bucket_ab'] == 9),
            'Category_ab'] = 'Predicted_Tractable__Medium_to_low_confidence'

        # self.out_df = self.out_df[(self.out_df['Top_bucket'] < 9 ) | (self.out_df['Top_bucket_ab'] < 10) ]

        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})


class Protac_buckets(object):
    '''
    Class for assigning genes to tractability buckets
    '''

    ##############################################################################################################
    #
    # Initial setup
    #
    #
    ##############################################################################################################

    def __init__(self, Pipeline_setup, database_url=None, ab_output=None):
        self.store_fetched = Pipeline_setup.store_fetched
        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # Cross referencing from Pipeline_setup, prevents repetition for antibody ot_tractability_pipeline

        self.id_xref = Pipeline_setup.id_xref

        # If antibody results are to be combined with small molecule results, append antibody columns to sm results
        # Otherwise, use the id_xref dataframe

        if ab_output is not None:
            go_data = self.id_xref[['ensembl_gene_id', 'go.CC']]
            self.out_df = ab_output.merge(go_data, how='outer', on='ensembl_gene_id')

        else:
            self.out_df = self.id_xref

        # ChEMBL currently not used

        # if database_url is None:
        #     database_url = os.getenv('CHEMBL_DB')
        #
        #
        # # Create ChEMBL DB connection
        # self.engine = create_engine(database_url)


    def _high_conf_locations(self, row):

        if isinstance(row['Uniprot_high_conf_loc'], str):
            if len(row['Uniprot_high_conf_loc']) == 0:
                row['Uniprot_high_conf_loc'] = []
            else: 
                row['Uniprot_high_conf_loc'] = eval(row['Uniprot_high_conf_loc'])

        if isinstance(row['GO_high_conf_loc'], str):
            if len(row['GO_high_conf_loc']) == 0:
                row['GO_high_conf_loc'] = []
            else: 
                row['GO_high_conf_loc'] = eval(row['GO_high_conf_loc'])

        try:
            len(row['Uniprot_high_conf_loc'])
        except TypeError:
            row['Uniprot_high_conf_loc'] = []

        try:
            len(row['GO_high_conf_loc'])
        except TypeError:
            row['GO_high_conf_loc'] = []



        if len(row['Uniprot_high_conf_loc']) == 0 and len(row['GO_high_conf_loc']) == 0 and row['PROTAC_location_Bucket'] == 5:
            return 5




        # locations = [x[0].lower().strip() for x in eval(row['Uniprot_high_conf_loc'])] + [x[0] for x in eval(row['GO_high_conf_loc'])]
        locations = [x[0].lower().strip() for x in row['Uniprot_high_conf_loc']] + [x[0] for x in row['GO_high_conf_loc']]
        accepted_locations = list(set(locations) & set(self.good_locations))
        grey_locations = list(set(locations) & set(self.grey_locations))
        # bad_locations = list(set(locations) & set(self.bad_locations))



        if len(accepted_locations) > 0:
            return 1

        elif len(grey_locations) > 0:
            return 3

        # elif len(bad_loactions) > 0:
        #     return 7


        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to
        # be bad

        elif row['PROTAC_location_Bucket'] == 6:
            return 7

        else:
            #print('locations',locations, 'accepted', accepted_locations, 'grey', grey_locations)

            return row['PROTAC_location_Bucket']




    def _med_conf_locations(self,row):

        if isinstance(row['Uniprot_med_conf_loc'], str):
            if len(row['Uniprot_med_conf_loc']) == 0:
                row['Uniprot_med_conf_loc'] = []
            else: 
                row['Uniprot_med_conf_loc'] = eval(row['Uniprot_med_conf_loc'])

        if isinstance(row['GO_med_conf_loc'], str):
            if len(row['GO_med_conf_loc']) == 0:
                row['GO_med_conf_loc'] = []
            else: 
                row['GO_med_conf_loc'] = eval(row['GO_med_conf_loc'])

        try:
            len(row['Uniprot_med_conf_loc'])
        except TypeError:
            row['Uniprot_med_conf_loc'] = []

        try:
            len(row['GO_med_conf_loc'])
        except TypeError:
            row['GO_med_conf_loc'] = []

        # should this be `if len(row['Uniprot_high_conf_loc']) == 0 and len(row['GO_high_conf_loc']) == 0 and row['PROTAC_location_Bucket'] == 5:`?
        if len(row['Uniprot_med_conf_loc']) == 0 and len(row['GO_med_conf_loc']) == 0:
            return 5



        #locations = [x[0].lower().strip() for x in eval(row['Uniprot_med_conf_loc'])] + [x[0] for x in eval(row['GO_med_conf_loc'])]
        locations = [x[0].lower().strip() for x in row['Uniprot_med_conf_loc']] + [x[0] for x in row['GO_med_conf_loc']]
        accepted_locations = list(set(locations) & set(self.good_locations))
        grey_locations = list(set(locations) & set(self.grey_locations))
        # bad_locations = list(set(locations) & set(self.bad_locations))


        if len(accepted_locations) > 0:
            return 2

        elif len(grey_locations) > 0:
            return 4

        # elif len(bad_loactions) > 0:
        #     return 6

        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to
        # be bad

        else:
            return 6



    def _PROTAC_location_bucket(self):

        ''''
        For PROTACs, only intracellular targets are suitable. Therefore, we will assign a score based on location to
        allow the targets to be filtered to those in the cytosol, nucleus or membrane with accessible portion.

        1 - High confidence good location
        2 - Med confidence good location
        3 - High confidence grey location
        4 - Med condifence grey location
        5 - Unknown location
        6 - Med confidence bad location
        7 - High confidence bad location
        '''

        self.good_locations = ['cytoplasm', 'cytosol', 'nucleus']
        self.grey_locations = ['membrane']
        # self.bad_locations = ['secreted']



        self.out_df['PROTAC_location_Bucket'] = 0

        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._med_conf_locations, axis=1)
        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._high_conf_locations, axis=1)



    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    # Clinical PROTAC targets
    #
    ##############################################################################################################

    def _assign_buckets_1_to_3(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        Currently, PROTACs are not labelled in ChEMBL, and at the time of writing, only Androgen Receptor
        (ENSG00000169083) has a phase 1 PROTAC. For now, this will be returned manually, but a decision should be made
        about how PROTACs are labelled in ChEMBL

        :return:
        '''

        self.out_df['Bucket_1_PROTAC'] = 0
        self.out_df['Bucket_2_PROTAC'] = 0
        self.out_df['Bucket_3_PROTAC'] = 0

        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000169083'), 'Bucket_3_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 4 and 5
    # Protein Turnover
    #
    ##############################################################################################################

    def _assign_bucket_4_and_5(self):
        '''
        Protein Turnover
        '''

        self.out_df['Bucket_4_PROTAC'] = 0
        self.out_df['Bucket_5_PROTAC'] = 0

        df = pd.read_csv(os.path.join(DATA_PATH, 'protein_half_life_hq.csv'))

        df = df.merge(self.out_df, right_on='symbol', left_on='gene_name', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max()
        df['Max_halflife'].fillna(-1, inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id', suffixes=['_drop', ''])

        self.out_df.loc[(self.out_df['Max_halflife'] >= 24), 'Bucket_4_PROTAC'] = 1
        self.out_df.loc[(self.out_df['Max_halflife'] > 10) & (self.out_df['Max_halflife'] < 24), 'Bucket_5_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 6
    # Known ubiquitination sites
    #
    ##############################################################################################################

    def _assign_bucket_6(self):
        '''
        Known ubiquitation sites
        '''

        ub_df = pd.read_csv(os.path.join(DATA_PATH, 'ubiquitination_sites.csv'))
        self.out_df = ub_df.merge(self.out_df, on='symbol', how='right')

        self.out_df['Bucket_6_PROTAC'] = 0
        self.out_df.loc[(self.out_df['number_of_ubiquitination_sites'] > 0), 'Bucket_6_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    # Predicted ubiquitination sites
    #
    ##############################################################################################################

    def _assign_bucket_7(self):
        '''
        Predicted ubiquitination sites

        '''
        self.out_df['Bucket_7_PROTAC'] = 0

    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Taregts mentioned in PROTAC literature
    #
    ##############################################################################################################

    def _search_papers(self):

        url = urllib2.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22proteolysis%20targeting%20chimera%22&resultType=lite&cursorMark=*&pageSize=1000&format=json")
        data = url.read()
        try: data = json.loads(data.decode())
        except UnicodeDecodeError: data = json.loads(data)
        df = pd.read_json(json.dumps(data['resultList']['result']), orient='records')

        return df[['authorString', 'id', 'issue',
                   'journalTitle', 'pmcid',
                   'pmid', 'pubType', 'pubYear', 'source', 'title', 'tmAccessionTypeList']]

    def _get_tagged_targets(self):

        articles = list(self.papers_df['search_id'].unique())

        # API only able to accept 8 IDs at a time
        n = 7
        chunks = ['&'.join(articles[i:i + n]) for i in range(0, len(articles), n)]
        df_lists = []
        tags_list = []
        for chunk in chunks:
            url_s = 'https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?{}&type=Gene_Proteins&format=JSON'.format(chunk)
            url = urllib2.urlopen(url_s)
            data = url.read()
            try: data = json.loads(data.decode())
            except UnicodeDecodeError: data = json.loads(data)
            annot_df = json_normalize(data,
                                      record_path='annotations')  # pd.read_json(json.dumps(data), orient='records')
            tags_df = json_normalize(data, record_path=['annotations', 'tags'])
            df_lists.append(annot_df)
            tags_list.append(tags_df)

            time.sleep(1.5)

        self.annotations = pd.concat(df_lists)
        self.tags = pd.concat(tags_list)
        self.annotations.reset_index(inplace=True)
        self.tags.reset_index(inplace=True)

        if self.store_fetched: 
            self.annotations.to_csv("{}/protac_pmc_annotations.csv".format(self.store_fetched), encoding='utf-8')
            self.tags.to_csv("{}/protac_pmc_tags.csv".format(self.store_fetched), encoding='utf-8')

    def _extract_uniprot(self, row):
        try:
            return row['uri'].split('/')[-1]
        except AttributeError:
            return row['uri']

    def _extract_id(self, row):
        try:
            short_id = row['id'].split('/')[-1].split('#')[0]
            return short_id
        except AttributeError:
            return None

    def _process_IDs(self):
        grouped_tags = self.tags.groupby('name').first()
        tagged_annotations = self.annotations.merge(grouped_tags, how='left', left_on='exact', right_on='name')
        tagged_annotations['accession'] = tagged_annotations.apply(self._extract_uniprot, axis=1)

        tagged_annotations['short_id'] = tagged_annotations.apply(self._extract_id, axis=1)

        # tagged_annotations['short_id']
        joined = tagged_annotations.merge(self.papers_df, left_on='short_id', right_on='id', how='inner')

        return joined[['accession', 'prefix', 'exact', 'postfix', 'section', 'full_id', 'journalTitle']]

    def _search_ID(self, row):
        return "articleIds={}%3A{}".format(row['source'], row['id'])

    def _full_ID(self, row):
        return "http://europepmc.org/abstract/{}/{}#eur...".format(row['source'], row['id'])

    def _assign_bucket_8(self):
        '''
        Mentioned in PROTAC literature
        '''

        self.papers_df = self._search_papers()
        
        if self.store_fetched: 
            self.papers_df.to_csv("{}/protac_pmc_papers.csv".format(self.store_fetched), encoding='utf-8')

        self.papers_df['search_id'] = self.papers_df.apply(self._search_ID, axis=1)
        self.papers_df['full_id'] = self.papers_df.apply(self._full_ID, axis=1)

        self._get_tagged_targets()

        tagged_targets_df = self._process_IDs()

        self.out_df = self.out_df.merge(tagged_targets_df, how='left', on='accession')

        self.out_df['Bucket_8_PROTAC'] = 0
        self.out_df.loc[(~self.out_df['full_id'].isna()), 'Bucket_8_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Small Molecule Tractable
    #
    ##############################################################################################################

    def _assign_bucket_9(self):
        '''
        Small molecule tractable
        '''

        self.out_df['Bucket_9_PROTAC'] = 0
        self.out_df.loc[(self.out_df['Top_bucket'] < 9), 'Bucket_9_PROTAC'] = 1

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    # def _clinical_precedence(self, s):
    #     return 1 * s['Bucket_1_ab'] + 0.7 * s['Bucket_2_ab'] + 0.2 * s['Bucket_3_ab']
    #
    # def _high_conf_pred(self, s):
    #     return 0.7 * s['Bucket_4_ab'] + 0.3 * s['Bucket_5_ab']
    #
    # def _med_conf_pred(self, s):
    #     return 0.4 * s['Bucket_6_ab'] + 0.25 * s['Bucket_7_ab'] + 0.25 * s['Bucket_8_ab'] + 0.1 * s['Bucket_9_ab']

    def _summarise_buckets(self):

        self.out_df['Top_bucket_PROTAC'] = 10
        for x in range(9, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_PROTAC'.format(x)] == 1), 'Top_bucket_PROTAC'] = x
            self.out_df['Bucket_{}_PROTAC'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_PROTAC'] = self.out_df['Bucket_1_PROTAC'] + self.out_df['Bucket_2_PROTAC'] + \
                                           self.out_df[
                                               'Bucket_3_PROTAC'] + self.out_df['Bucket_4_PROTAC'] + self.out_df[
                                               'Bucket_5_PROTAC'] + self.out_df['Bucket_6_PROTAC'
                                           ] + self.out_df['Bucket_7_PROTAC'] + self.out_df['Bucket_8_PROTAC'] + \
                                           self.out_df[
                                               'Bucket_9_PROTAC']

        self.out_df.set_index('ensembl_gene_id')

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._PROTAC_location_bucket()
        self._assign_buckets_1_to_3()
        self._assign_bucket_4_and_5()
        self._assign_bucket_6()
        self._assign_bucket_7()
        self._assign_bucket_8()
        self._assign_bucket_9()

        self._summarise_buckets()

        # try:

        # self.out_df.index = self.out_df['ensembl_gene_id']
        self.out_df = self.out_df.groupby('ensembl_gene_id').first()

        # Columns to keep. This includes columns from the small molecule pipeline
        self.out_df = self.out_df[['accession', 'symbol',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4',
                                   'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_sum', 'Top_bucket', 'Category',
                                   'Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable',
                                   'PDB_Known_Ligand',
                                   'ensemble', 'High_Quality_ChEMBL_compounds',
                                   'Small_Molecule_Druggable_Genome_Member',
                                   'Bucket_1_ab', 'Bucket_2_ab', 'Bucket_3_ab', 'Bucket_4_ab',
                                   'Bucket_5_ab', 'Bucket_6_ab', 'Bucket_7_ab',
                                   'Bucket_8_ab', 'Bucket_9_ab', 'Bucket_sum_ab', 'Top_bucket_ab',
                                   'Clinical_Precedence_ab', 'Predicted_Tractable__High_confidence', 'Predicted_Tractable__Medium_to_low_confidence', 'Category_ab',
                                   'Uniprot_high_conf_loc', 'GO_high_conf_loc',
                                   'Uniprot_med_conf_loc',
                                   'GO_med_conf_loc', 'Transmembrane', 'Signal_peptide', 'HPA_main_location',
                                   'Bucket_1_PROTAC', 'Bucket_2_PROTAC', 'Bucket_3_PROTAC', 'Bucket_4_PROTAC',
                                   'Bucket_5_PROTAC', 'Bucket_6_PROTAC', 'Bucket_7_PROTAC',
                                   'Bucket_8_PROTAC', 'Bucket_9_PROTAC', 'Bucket_sum_PROTAC', 'Top_bucket_PROTAC',
                                   'Bcell_mean', 'NKcell_mean', 'Hepatocytes_mean', 'MouseNeuorons_mean',
                                   'Max_halflife',
                                   'number_of_ubiquitination_sites',
                                   'full_id', 'PROTAC_location_Bucket'
                                   ]]
        self.out_df.sort_values(['Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable'],
                                ascending=[False, False, False], inplace=True)


        self.out_df = self.out_df[(self.out_df['Top_bucket'] < 9) | (self.out_df['Top_bucket_ab'] < 10) | (
                    self.out_df['Top_bucket_PROTAC'] < 10)]

        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})
