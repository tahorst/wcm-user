'''
Replaces expression in wcm data with sequencing from Gene Wei Li.
Can use seal or rsem analysis from wcm by setting WCM_SEQ below.

Output:
	tsv with updated expression saved in out/
'''

from __future__ import division
import csv
import os

import numpy as np

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli


SCRIPT_DIR = os.path.dirname(__file__)
GENE_WEI_LI_FILE = os.path.join(SCRIPT_DIR, 'rna_seq_gene_wei_li.tsv')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'out')
if not os.path.exists(OUTPUT_DIR):
	os.mkdir(OUTPUT_DIR)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'updated_seq.tsv')

WCM_SEQ = 'rnaseq_rsem_tpm_mean'
# WCM_SEQ = 'rnaseq_seal_rpkm_mean'  # Uncomment for expression from seal analysis


if __name__ == '__main__':
	# Load data
	raw_data = KnowledgeBaseEcoli()

	wcm_seq = {seq['Gene']: seq['M9 Glucose minus AAs'] for seq in getattr(raw_data.rna_seq_data, WCM_SEQ)}

	with open(GENE_WEI_LI_FILE) as f:
		reader = csv.reader(f, delimiter='\t')
		gwl = {line[0]: float(line[1]) for line in reader}

	gwl_to_wcm = {gene['symbol']: gene['id'] for gene in raw_data.genes}

	# Get counts of Gene Wei Li genes in both datasets to scale to the wcm expression
	valid_gwl_genes = [gene for gene in gwl if gene in gwl_to_wcm]
	wcm_counts = np.sum([wcm_seq[gwl_to_wcm[gene]] for gene in valid_gwl_genes])
	valid_gwl_counts = np.sum([gwl[gene] for gene in valid_gwl_genes])
	rescale = wcm_counts / valid_gwl_counts

	# Replace genes with Gene Wei Li expression
	updated_counts = {}
	updated_counts.update(wcm_seq)
	updated_counts.update({gwl_to_wcm[gene]: gwl[gene] * rescale for gene in valid_gwl_genes})

	# Normalize final counts to transcripts per million
	tpm_scale = 1000000 / np.sum(updated_counts.values())

	# Save updated expression
	with open(OUTPUT_FILE, 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		for gene, exp in updated_counts.items():
			writer.writerow([gene, exp * tpm_scale])
