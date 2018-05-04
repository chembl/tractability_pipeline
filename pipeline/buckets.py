import json
import sys
import time
import zipfile
import io

import mygene
import pandas as pd
from sqlalchemy import create_engine

from pipeline.sm_queries import *
from pipeline.ab_queries import *

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
    from functools import reduce
else:
    import urllib2


class Pipeline_setup(object):
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
                              as_dataframe=True, fields='uniprot.Swiss-Prot,entrezgene,pdb,pfam,go,interpro',
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

    def __init__(self, Pipeline_setup, database_url, ligand_filter=list()):

        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # Cross referencing from Pipeline_setup, prevents repetition for antibody pipeline

        self.id_xref = Pipeline_setup.id_xref

        # Load lists for PDB ligand filters (i.e. to remove sugars, solvents etc)

        if len(ligand_filter) == 0:
            # Default ligand filter
            # Load unwanted PDB ligands list from text files

            with urllib2.urlopen(
                    'https://gist.githubusercontent.com/Cradoux/8d16493e19de89f7134d0a800e2dd06c/raw/5c697de2ebdbce6694df875e5e8338dbfed68171/organic_solvents.txt') as solvent_file:
                solvents = [a.decode("utf-8").split('\t')[0] for a in solvent_file]

            with urllib2.urlopen(
                    'https://gist.githubusercontent.com/Cradoux/adb2e77918828685d3b60418dc3ad7ca/raw/2a6a0b39bb1b1701536fa8abbc6b728615b2aa6a/sugars.txt') as sugar_file:
                sugars = [a.decode("utf-8").split('\t')[0] for a in sugar_file]

            with urllib2.urlopen(
                    'https://gist.githubusercontent.com/Cradoux/3bd1fbf834eea6d9ee2239fcb139cade/raw/d7fd1fd2a655abcc8e917a8e4a03a97ed1aa2b2f/cofactors.txt') as cofactor_file:
                cofactors = [a.decode("utf-8").split('\t')[0] for a in cofactor_file]

            ligand_filter = solvents + sugars + cofactors

            self.ligand_filter = ligand_filter

        # Dataframes for buckets

        # All chembl data loaded into here
        self.all_chembl_targets = None

        # Unique list of PDB codes:
        self.pdb_list = []

        # Map back to Uniprot accession
        self.pdb_map = {}

        # Create ChEMBL DB connection
        self.engine = create_engine(database_url)

        # URL for PDBe web services

        self.SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

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
            from CHEMBL_23.target_dictionary td, CHEMBL_23.binding_sites bs
            where td.tid = bs.tid and td.tid IN {}
            
            '''.format(tuple(chunk))
            df_list.append(pd.read_sql_query(q, self.engine))

        # Merge will set those with unknown binding site as NAN
        binding_site_info = pd.concat(df_list)
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
            from CHEMBL_23.component_sequences cs, CHEMBL_23.site_components sc
            where cs.component_id = sc.component_id
            and cs.accession in {}'''.format(tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2)
        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc])

    def _assign_buckets_1_to_3(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        Possible duplication of effort of the OT known_drug score

        :return:
        '''

        def other_func(x):
            return tuple(x)

        self._search_chembl_clinical()
        self._process_protein_complexes()

        self.gene_xref = self.id_xref[['accession','ensembl_gene_id']]

        self.out_df = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')

        self.out_df.drop(['component_id', 'drug_name', 'ref_id', 'ref_type', 'tid', 'molregno',
                          'parent_molregno', 'ref_url'], axis=1, inplace=True)

        f = {x: 'first' for x in self.out_df.columns}
        f['max_phase'] = 'max'
        f['pref_name'] = other_func
        f['moa_chembl'] = other_func

        self.out_df = self.out_df.groupby(['ensembl_gene_id']).agg(f)

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
        full_url = "%s/%s/?pretty=%s" % (self.SERVER_URL, url, str(pretty).lower())

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
        df = pd.read_csv('drugebility_scores.csv')

        df = df.merge(self.gene_xref, on='accession', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max()
        df['ensemble'].fillna(-1, inplace=True)

        self.out_df = df.merge(self.out_df, on='ensembl_gene_id', how='right')
        self.out_df['Bucket_5'] = 0
        self.out_df['Bucket_6'] = 0



        self.out_df.loc[(self.out_df['ensemble'] >= 0.7),'Bucket_5'] = 1
        self.out_df.loc[(self.out_df['ensemble'] > 0) & (self.out_df['ensemble'] < 0.7),'Bucket_6'] = 1

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
                                     pd.read_sql_query(Emax_q, self.engine)])

        # For faster testing
        # self.activities = pd.read_sql_query(pchembl_q, self.engine)

    def _structural_alerts(self):
        q = '''
        select distinct cs.canonical_smiles, csa.alert_id /*, sa.alert_name, sas.set_name */
        from CHEMBL_23.compound_structures cs,
            CHEMBL_23.compound_structural_alerts csa,
            CHEMBL_23.structural_alerts sa,
            CHEMBL_23.molecule_dictionary md,
            CHEMBL_23.structural_alert_sets sas
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


        df = self.activities[(self.activities['pfi'] <= 7) & (self.activities['alert_id'] <= 2)]

        f = {x: 'first' for x in df.columns}
        f['canonical_smiles'] = 'count'
        print(df)
        df2 = df.groupby('accession').agg(f).reset_index(drop=True)
        #df2['target_chembl_id'] = df.groupby('accession', as_index=False)['target_chembl_id'].first()
        df2 = df2[['accession','canonical_smiles','target_chembl_id']]
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
        df = pd.read_csv('druggable_genome.csv')
        df = df[['ensembl_gene_id','small_mol_druggable']]
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

    def _regroup_buckets(self):

        bucket_li = [self.b1, self.b2, self.b3, self.b4, self.b5, self.b6, self.b7, self.b8, self.b9]
        i = 0
        for b in bucket_li:
            i += 1
            try:
                b['Bucket'] = i
                b['Bucket_{}'.format(i)] = 1
                b.dropna(axis=1,how='all', inplace=True)
            except TypeError:
                # In case of empty bucket
                continue

        i = 0
        for b in bucket_li:
            i += 1
            try:
                b['Bucket_{}'.format(i)].fillna(0, inplace=True)
            except TypeError:
                continue
            except KeyError:
                continue



        df = reduce(lambda left, right: pd.merge(left, right, on='ensembl_gene_id', how='outer'),
                    [b for b in bucket_li if b is not None])

        # df = df.drop(['Unnamed: 0.1', '_id', '_score', 'acd_logd', 'adme_gene', 'alert_id', 'aromatic_rings',
        #               'canonical_smiles', 'chr_b37', 'compound_chembl_id', 'description',
        #               'end_b37', 'molregno', 'no_of_gwas_regions', 'parent_molregno', 'pfi', 'px_number',
        #               'qed_weighted', 'start_b37', 'strand', 'year'], axis = 1)

        # df['Bucket_sum'] = df['Bucket_1'] + df['Bucket_2'] + df['Bucket_3'] + df['Bucket_4'] + df['Bucket_5'] + df[
        #     'Bucket_6'] + df['Bucket_7'] + df['Bucket_8']

        return df

    def _summarise_buckets(self):

        self.out_df['Bucket_sum'] = self.out_df['Bucket_1'] + self.out_df['Bucket_2']+ self.out_df[
            'Bucket_3']+ self.out_df['Bucket_4'] + self.out_df['Bucket_5'] + self.out_df['Bucket_6'] + self.out_df[
            'Bucket_7'] + self.out_df['Bucket_8']

        self.out_df['Top_bucket'] = 10
        self.out_df.loc[(self.out_df['Bucket_8'] == 1), 'Top_bucket'] = 8
        self.out_df.loc[(self.out_df['Bucket_7'] == 1), 'Top_bucket'] = 7
        self.out_df.loc[(self.out_df['Bucket_6'] == 1), 'Top_bucket'] = 6
        self.out_df.loc[(self.out_df['Bucket_5'] == 1), 'Top_bucket'] = 5
        self.out_df.loc[(self.out_df['Bucket_4'] == 1), 'Top_bucket'] = 4
        self.out_df.loc[(self.out_df['Bucket_3'] == 1), 'Top_bucket'] = 3
        self.out_df.loc[(self.out_df['Bucket_2'] == 1), 'Top_bucket'] = 2
        self.out_df.loc[(self.out_df['Bucket_1'] == 1), 'Top_bucket'] = 1


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
        #self._assign_bucket_9()
        self._summarise_buckets()

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

    def __init__(self, Pipeline_setup, database_url, sm_output = None):

        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # Cross referencing from Pipeline_setup, prevents repetition for antibody pipeline

        self.id_xref = Pipeline_setup.id_xref

        # If antibody results are to be combined with small molecule results, append antibody columns to sm results
        # Otherwise, use the id_xref dataframe

        if sm_output is not None:
            self.out_df = sm_output
        else:
            self.out_df = self.id_xref

        # Load lists for PDB ligand filters (i.e. to remove sugars, solvents etc)

        # Dataframes for buckets

        # All chembl data loaded into here
        self.all_chembl_targets = None

        # Unique list of PDB codes:
        self.pdb_list = []

        # Map back to Uniprot accession
        self.pdb_map = {}

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
            from CHEMBL_23.component_sequences cs,
                CHEMBL_23.site_components sc
            where cs.component_id = sc.component_id
            and sc.site_id in {}'''.format(tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2)
        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc])

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

        df = self.all_chembl_targets.merge(self.id_xref, how='outer', on='accession')
        df = df.groupby('ensembl_gene_id', as_index=False).first()
        df['max_phase'].fillna(0, inplace=True)

        # Bucket 1: Targets with Phase 4
        self.b1 = df[df['max_phase'] == 4]

        # Bucket 2: Targets >= Phase 2
        self.b2 = df[(df['max_phase'] < 4) & (df['max_phase'] >= 2)]

        # Bucket 3: Targets with lead op or preclinical SM ################## Assume phase 0-2 for now
        self.b3 = df[(df['max_phase'] < 2) & (df['max_phase'] > 0)]

    ##############################################################################################################
    #
    # Functions relating to buckets 4 to
    #
    #
    ##############################################################################################################

    def _assign_bucket_4(self):
        '''
        Uniprot (loc): Targets in "Cell membrane" or "Secreted", high confidence
        '''

        print(self.id_xref)

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
        pass

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    #
    #
    ##############################################################################################################

    def _assign_bucket_7(self):
        '''
        Does the target have ligands in ChEMBL (PFI <=7, SMART hits <= 2, scaffolds >= 2)
        '''
        pass

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
        pass

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
        Search to see whether targets have chemical patents in last 5 years

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
                               "Cell cycle dependency", "GO id"])

        self.b9 = reliable[reliable['main_location'].str.contains("Plasma membrane", na=False)]

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _regroup_buckets(self):

        bucket_li = [self.b1, self.b2, self.b3, self.b4, self.b5, self.b6, self.b7, self.b8, self.b9]
        i = 0
        for b in bucket_li:
            i += 1
            try:
                b['Bucket'] = i
                b['Bucket_{}'.format(i)] = 1
            except TypeError:
                # In case of empty bucket
                continue

        out_df = pd.concat(bucket_li)

        i = 0
        for b in bucket_li:
            i += 1
            try:
                b['Bucket_{}'.format(i)].fillna(0, inplace=True)
            except TypeError:
                continue

        return out_df

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
        self._assign_bucket_9()

        return self._regroup_buckets()
