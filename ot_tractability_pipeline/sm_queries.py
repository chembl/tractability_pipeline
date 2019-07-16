# The length of the SQL queries makes the Python code far too long, moved to separate file

import os

try:
    chembl_number = int(os.getenv('CHEMBL_VERSION'))
    CHEMBL_VERSION = "CHEMBL_{}".format(chembl_number)
except TypeError:
    print('''
    Please set the CHEMBL version '''
          )
    raise

            
chembl_clinical_small_mol = """
        select distinct mh.parent_molregno, 
            md.molregno, 
            md.pref_name
        from {0}.molecule_dictionary md,
            {0}.molecule_hierarchy mh
        where md.molregno = mh.molregno
        
        and md.molecule_type = 'Small molecule'
            """.format(CHEMBL_VERSION)

chembl_clinical_targets = """
        select distinct mh.parent_molregno, 
            md.pref_name as drug_name, 
            dt.development_phase as max_phase,
            td.tid,
            td.chembl_id as target_chembl_id,
            td.pref_name as target_name, 
            td.target_type,
            cs.accession as accession,
            dm.action_type as moa_chembl,
            mr.ref_type,
            mr.ref_id,
            mr.ref_url
        from {0}.molecule_dictionary md,
            {0}.molecule_hierarchy mh,
            {0}.drug_mechanism dm,
            {0}.target_dictionary td,
            {0}.target_components tc,
            {0}.component_sequences cs,
            {0}.mechanism_refs mr,
            CHEMBL_APP.drug_targets_mv dt
        where md.molregno = dm.molregno
        and md.molregno = mh.molregno
        and dm.tid = td.tid
        and td.tid = tc.tid
        and dt.tid = tc.tid
        and tc.component_id = cs.component_id
        and dm.mec_id = mr.mec_id
        
        and md.molecule_type = 'Small molecule'
        and td.tax_id = 9606
        and td.target_type like '%PROTEIN%'""".format(CHEMBL_VERSION)

pchembl_q = '''
        select distinct td.chembl_id as target_chembl_id, 
            td.tid,
            cs.accession, 
            pf.protein_class_desc, 
            md.molregno, 
            md.chembl_id as compound_chembl_id,
            cst.canonical_smiles,
            cp.acd_logd,
            cp.aromatic_rings,
            cp.qed_weighted,
            d.year
        from {0}.target_dictionary td,
            {0}.target_components tc,
            {0}.component_sequences cs,
            {0}.component_class cc,
            {0}.protein_classification pf,
            {0}.assays ass,
            {0}.activities act,
            {0}.molecule_dictionary md,
            {0}.compound_structures cst,
            {0}.compound_properties cp,
            {0}.docs d
        where td.tid = tc.tid
        and tc.component_id = cs.component_id
        and cs.component_id = cc.component_id
        and cc.protein_class_id = pf.protein_class_id
        and td.tid = ass.tid
        and ass.assay_id = act.assay_id
        and act.molregno = md.molregno
        and md.molregno = cst.molregno
        and md.molregno = cp.molregno
        and act.doc_id = d.doc_id
        and td.tax_id = 9606
        and td.target_type like '%PROTEIN%'
        and ass.assay_type = 'B'
        and ass.relationship_type = 'D'
        and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
        and act.potential_duplicate = 0
        and (act.pchembl_value is not NULL and pchembl_value >= 5.5)
        and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
        and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))
        '''.format(CHEMBL_VERSION)
nm_q = '''
select distinct td.chembl_id as target_chembl_id, 
    td.tid, 
    cs.accession, 
    pf.protein_class_desc, 
    md.molregno, 
    md.chembl_id as compound_chembl_id,
    cst.canonical_smiles,
    cp.acd_logd,
    cp.aromatic_rings,
    cp.qed_weighted,
    d.year
from {0}.target_dictionary td,
    {0}.target_components tc,
    {0}.component_sequences cs,
    {0}.component_class cc,
    {0}.protein_classification pf,
    {0}.assays ass,
    {0}.activities act,
    {0}.molecule_dictionary md,
    {0}.compound_structures cst,
    {0}.compound_properties cp,
    {0}.docs d
where td.tid = tc.tid
and tc.component_id = cs.component_id
and cs.component_id = cc.component_id
and cc.protein_class_id = pf.protein_class_id
and td.tid = ass.tid
and ass.assay_id = act.assay_id
and act.molregno = md.molregno
and md.molregno = cst.molregno
and md.molregno = cp.molregno
and act.doc_id = d.doc_id
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
and ass.assay_type = 'B'
and ass.relationship_type = 'D'
and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
and act.potential_duplicate = 0
and (act.pchembl_value is NULL and act.standard_type in ('AC50', 'CC50', 'EC50', 'GI50', 'IC50', 'IC90', 'IC95', 'IC99', 'Kd', 'Ki', 'LC50', 'MIC', 'MIC50', 'Potency', 'Kinact', 'KB', 'Activity', 'Ke', 'KA', 'IC100') and act.standard_units = 'nM' and act.standard_value <= 3000 and act.standard_relation in ('<', '='))
and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))'''.format(CHEMBL_VERSION)

km_kon_q = '''
select distinct td.chembl_id as target_chembl_id, 
    td.tid, 
    cs.accession, 
    pf.protein_class_desc, 
    md.molregno, 
    md.chembl_id as compound_chembl_id,
    cst.canonical_smiles,
    cp.acd_logd,
    cp.aromatic_rings,
    cp.qed_weighted,
    d.year
from {0}.target_dictionary td,
    {0}.target_components tc,
    {0}.component_sequences cs,
    {0}.component_class cc,
    {0}.protein_classification pf,
    {0}.assays ass,
    {0}.activities act,
    {0}.molecule_dictionary md,
    {0}.compound_structures cst,
    {0}.compound_properties cp,
    {0}.docs d
where td.tid = tc.tid
and tc.component_id = cs.component_id
and cs.component_id = cc.component_id
and cc.protein_class_id = pf.protein_class_id
and td.tid = ass.tid
and ass.assay_id = act.assay_id
and act.molregno = md.molregno
and md.molregno = cst.molregno
and md.molregno = cp.molregno
and act.doc_id = d.doc_id
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
and ass.assay_type = 'B'
and ass.relationship_type = 'D'
and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
and act.potential_duplicate = 0
and (act.pchembl_value is NULL and act.standard_type in ('k_on', 'Km'))
and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))
'''.format(CHEMBL_VERSION)

D_Tm_q = '''
select distinct td.chembl_id as target_chembl_id, 
    td.tid, 
    cs.accession, 
    pf.protein_class_desc, 
    md.molregno, 
    md.chembl_id as compound_chembl_id,
    cst.canonical_smiles,
    cp.acd_logd,
    cp.aromatic_rings,
    cp.qed_weighted,
    d.year
from {0}.target_dictionary td,
    {0}.target_components tc,
    {0}.component_sequences cs,
    {0}.component_class cc,
    {0}.protein_classification pf,
    {0}.assays ass,
    {0}.activities act,
    {0}.molecule_dictionary md,
    {0}.compound_structures cst,
    {0}.compound_properties cp,
    {0}.docs d
where td.tid = tc.tid
and tc.component_id = cs.component_id
and cs.component_id = cc.component_id
and cc.protein_class_id = pf.protein_class_id
and td.tid = ass.tid
and ass.assay_id = act.assay_id
and act.molregno = md.molregno
and md.molregno = cst.molregno
and md.molregno = cp.molregno
and act.doc_id = d.doc_id
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
and ass.assay_type = 'B'
and ass.relationship_type = 'D'
and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
and act.potential_duplicate = 0
and (act.pchembl_value is NULL and act.standard_type in ('deltaTm', 'Tm') and act.standard_units = 'degrees C' and act.standard_value >= 2 and act.standard_relation in ('>', '='))
and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))
'''.format(CHEMBL_VERSION)

residual_act_q = '''
select distinct td.chembl_id as target_chembl_id, 
    td.tid, 
    cs.accession, 
    pf.protein_class_desc, 
    md.molregno, 
    md.chembl_id as compound_chembl_id,
    cst.canonical_smiles,
    cp.acd_logd,
    cp.aromatic_rings,
    cp.qed_weighted,
    d.year
from {0}.target_dictionary td,
    {0}.target_components tc,
    {0}.component_sequences cs,
    {0}.component_class cc,
    {0}.protein_classification pf,
    {0}.assays ass,
    {0}.activities act,
    {0}.molecule_dictionary md,
    {0}.compound_structures cst,
    {0}.compound_properties cp,
    {0}.docs d
where td.tid = tc.tid
and tc.component_id = cs.component_id
and cs.component_id = cc.component_id
and cc.protein_class_id = pf.protein_class_id
and td.tid = ass.tid
and ass.assay_id = act.assay_id
and act.molregno = md.molregno
and md.molregno = cst.molregno
and md.molregno = cp.molregno
and act.doc_id = d.doc_id
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
and ass.assay_type = 'B'
and ass.relationship_type = 'D'
and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
and act.potential_duplicate = 0
and (act.pchembl_value is NULL and act.standard_type = 'Residual activity' and act.standard_units = '%' and act.standard_value <= 10 and act.standard_relation in ('<', '='))
and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))
'''.format(CHEMBL_VERSION)

Imax_q = '''
select distinct td.chembl_id as target_chembl_id, 
    td.tid, 
    cs.accession, 
    pf.protein_class_desc, 
    md.molregno, 
    md.chembl_id as compound_chembl_id,
    cst.canonical_smiles,
    cp.acd_logd,
    cp.aromatic_rings,
    cp.qed_weighted,
    d.year
from {0}.target_dictionary td,
    {0}.target_components tc,
    {0}.component_sequences cs,
    {0}.component_class cc,
    {0}.protein_classification pf,
    {0}.assays ass,
    {0}.activities act,
    {0}.molecule_dictionary md,
    {0}.compound_structures cst,
    {0}.compound_properties cp,
    {0}.docs d
where td.tid = tc.tid
and tc.component_id = cs.component_id
and cs.component_id = cc.component_id
and cc.protein_class_id = pf.protein_class_id
and td.tid = ass.tid
and ass.assay_id = act.assay_id
and act.molregno = md.molregno
and md.molregno = cst.molregno
and md.molregno = cp.molregno
and act.doc_id = d.doc_id
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
and ass.assay_type = 'B'
and ass.relationship_type = 'D'
and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
and act.potential_duplicate = 0
and (act.pchembl_value is NULL and act.standard_type in ('Activity', 'Imax') and act.standard_units = '%' and act.standard_value >= 70 and act.standard_relation in ('>', '='))
and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))
'''.format(CHEMBL_VERSION)

Emax_q = '''
select distinct td.chembl_id as target_chembl_id, 
    td.tid, 
    cs.accession, 
    pf.protein_class_desc, 
    md.molregno, 
    md.chembl_id as compound_chembl_id,
    cst.canonical_smiles,
    cp.acd_logd,
    cp.aromatic_rings,
    cp.qed_weighted,
    d.year
from {0}.target_dictionary td,
    {0}.target_components tc,
    {0}.component_sequences cs,
    {0}.component_class cc,
    {0}.protein_classification pf,
    {0}.assays ass,
    {0}.activities act,
    {0}.molecule_dictionary md,
    {0}.compound_structures cst,
    {0}.compound_properties cp,
    {0}.docs d
where td.tid = tc.tid
and tc.component_id = cs.component_id
and cs.component_id = cc.component_id
and cc.protein_class_id = pf.protein_class_id
and td.tid = ass.tid
and ass.assay_id = act.assay_id
and act.molregno = md.molregno
and md.molregno = cst.molregno
and md.molregno = cp.molregno
and act.doc_id = d.doc_id
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
and ass.assay_type = 'B'
and ass.relationship_type = 'D'
and (act.data_validity_comment is NULL or act.data_validity_comment = 'Manually validated')
and act.potential_duplicate = 0
and (act.pchembl_value is NULL and act.standard_type in ('Emax', 'Efficacy') and act.standard_units = '%' and act.standard_value >= 120 and act.standard_relation in ('>', '='))
and md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
and ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') and cp.mw_freebase <= 1500) or (md.molecule_type = 'Small molecule'))
'''.format(CHEMBL_VERSION)
