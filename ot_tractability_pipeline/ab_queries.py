from ot_tractability_pipeline.sm_queries import CHEMBL_VERSION

chembl_clinical_ab = '''
select distinct mh.parent_molregno, di.efo_id, di.efo_term, di.max_phase_for_ind
from {0}.molecule_dictionary md,
	{0}.molecule_hierarchy mh,
	{0}.drug_indication di
where md.molregno = mh.molregno
and md.molregno = di.molregno

and md.molecule_type = 'Antibody'
'''.format(CHEMBL_VERSION)

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
from {0}.molecule_dictionary md,
	{0}.molecule_hierarchy mh,
	{0}.drug_mechanism dm,
	{0}.target_dictionary td,
	{0}.target_components tc,
	{0}.component_sequences cs,
	{0}.mechanism_refs mr
where md.molregno = dm.molregno
and md.molregno = mh.molregno
and dm.tid = td.tid
and td.tid = tc.tid
and tc.component_id = cs.component_id
and dm.mec_id = mr.mec_id

and md.molecule_type = 'Antibody'
and td.tax_id = 9606
and td.target_type like '%PROTEIN%'
'''.format(CHEMBL_VERSION)