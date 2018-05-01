from pipeline.buckets import Small_molecule_buckets
import pandas as pd

t = Small_molecule_buckets()
t.parse_uniprot_xml('uniprot_sprot.xml')
