from pipeline.buckets import *
import csv


class Accession_Picker(object):

    def __init__(self, ensembl_list):
        self.need_checking = []
        self.mg = mygene.MyGeneInfo()
        self.exp_codes = ['EXP','IDA','IPI','IMP','IGI','IEP']
        self.htp_codes = ['HTP','HDA','HMP','HGI','HEP']
        self.computational_analysis_codes = ['ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA']
        self.author_statement_codes = ['TAS','NAS']
        self.curator_statement_codes = ['IC','ND']
        self.electionic_annotation_codes = ['IEA']
        self.ensembl_list = ensembl_list

    def best_evidence(self):

        results = self.mg.querymany(list(set(self.need_checking)), scopes='uniprot',
                              as_dataframe=True, fields='ensembl,entrezgene,pdb,pfam,go,interpro',
                              species='human', returnall=True)
        print(results)

        results['out'].to_csv('uniprot_evidence.csv', quoting=csv.QUOTE_NONNUMERIC)


    def _uniprot_primary_only(self, a):
        '''
        If multiple uniprot IDs, only take primary (assume first in list is primary) <== NEEDS CHECKING
        :return:
        '''
        if isinstance(a, dict):
            s = a['Swiss-Prot']

            if isinstance(s, list):

                self.need_checking += s
            return s
        else:
            return a

    def get_xref(self):

        results = self.mg.getgenes(list(self.ensembl_list), scopes='ensembl',
                              as_dataframe=True, fields='uniprot.Swiss-Prot,entrezgene,pdb,pfam,go,interpro',
                              species='human', returnall=True)

        id_xref = results
        results['uniprot'] = id_xref['uniprot'].apply(self._uniprot_primary_only)
        self.best_evidence()


ot_data = pd.read_csv('all_targets.csv', encoding='utf-8')
ensembl_id_list = list(ot_data['ensembl_gene_id'].unique())

ap = Accession_Picker(ensembl_id_list)
ap.get_xref()