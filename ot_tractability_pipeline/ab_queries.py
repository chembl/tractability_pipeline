chembl_clinical_ab = '''
select distinct mh.parent_molregno, di.efo_id, di.efo_term, di.max_phase_for_ind
from CHEMBL_23.molecule_dictionary md,
	CHEMBL_23.molecule_hierarchy mh,
	CHEMBL_23.drug_indication di
where md.molregno = mh.molregno
and md.molregno = di.molregno
and md.therapeutic_flag = 1
and md.molecule_type = 'Antibody'
'''

chembl_clinical_ab_targets = '''
select distinct mh.parent_molregno, 
	md.pref_name as drug_name, 
	md.max_phase,
	td.tid,
	td.pref_name as target_name, 
	td.target_type,
	cs.accession as accession,
	dm.action_type as moa_chembl,
	mr.ref_type,
	mr.ref_id,
	mr.ref_url,
	dm.site_id
from CHEMBL_23.molecule_dictionary md,
	CHEMBL_23.molecule_hierarchy mh,
	CHEMBL_23.drug_mechanism dm,
	CHEMBL_23.target_dictionary td,
	CHEMBL_23.target_components tc,
	CHEMBL_23.component_sequences cs,
	CHEMBL_23.mechanism_refs mr
where md.molregno = dm.molregno
and md.molregno = mh.molregno
and dm.tid = td.tid
and td.tid = tc.tid
and tc.component_id = cs.component_id
and dm.mec_id = mr.mec_id
and md.therapeutic_flag = 1
and md.molecule_type = 'Antibody'
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
'''