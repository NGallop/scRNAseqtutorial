import requests
import json
import pandas as pd
import scanpy as sc

sc.settings.figdir = "plots/data_investigation/"

response=requests.get("https://www.ebi.ac.uk/gwas/rest/api/studies/GCST003588/associations")
gwas_json = response.json()

associations_data = response.json()['_embedded']['associations']
gwas = pd.json_normalize(associations_data, 'loci')
risk_alleles = gwas.explode('strongestRiskAlleles').explode('authorReportedGenes').reset_index(drop=True)

risk_alleles_df = pd.json_normalize(risk_alleles['strongestRiskAlleles'])
reportedgenes_df = pd.json_normalize(risk_alleles['authorReportedGenes'])

gwasgenes=pd.concat([risk_alleles_df, reportedgenes_df], axis=1)
gwasgenes = gwasgenes[['riskAlleleName', 'geneName']]

genes=[gene for gene in gwasgenes['geneName'].tolist() if gene !='Intergenic']
