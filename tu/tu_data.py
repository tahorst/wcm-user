'''
Uses transcription unit information to compare deg rates and expression within operons.
Also produces counts for how many genes and TUs are in each operon.
Should run with the wcEcoli repo in PYTHONPATH

Outputs multiple files:
	tu_deg.tsv - deg rate for each gene within a TU, grouped by operon
	tu_exp.tsv - exp for each gene within a TU, grouped by operon
	gene_to_tu_count.tsv - gene and TU counts per operon
'''

import numpy as np
import json

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.simulation_data import SimulationDataEcoli

raw_data = KnowledgeBaseEcoli()
sim_data = SimulationDataEcoli()
sim_data.initialize(raw_data)

transcription = sim_data.process.transcription

deg = {rna['id']: rna['degRate'] for rna in transcription.rnaData}
exp = {id: e for id, e in zip(transcription.rnaData['id'], transcription.rnaExpression['basal'])}

genes_included = {id: False for id in exp}

sets = []

# loop over each gene to get all operons (TUs that share common genes)
for gene in genes_included:
	# skip genes that have already been counted
	if genes_included[gene]:
		continue

	genes = set([gene])
	length = 0

	# continue until no more genes are added from shared TUs
	while len(genes) != length:
		length = len(genes)
		tus = np.hstack(transcription.getTranscriptionUnitWithRna(list(genes)))
		genes = set(np.hstack(transcription.getRnaInTranscriptionUnit(tus)))

	# mark gene to prevent future counting
	for g in genes:
		genes_included[g] = True

	# save operons for file output
	sets.append(genes)

# save deg rates for all genes in a TU, grouped by operons
with open('tu_deg.tsv', 'w') as f:
	for genes in sets:
		tus = set(np.hstack(transcription.getTranscriptionUnitWithRna(list(genes))))
		singles = np.sum([1 if len(tu) == 1 else 0 for tu in transcription.getRnaInTranscriptionUnit(list(tus))])
		f.write('{} genes in {} TUs, with {} single genes\n'.format(len(genes), len(tus), singles))
		for tu in transcription.getRnaInTranscriptionUnit(list(tus)):
			f.write(json.dumps(tu.tolist()))
			for gene in tu:
				f.write('\t{}'.format(deg[gene]))
			f.write('\n')

# save expression for all genes in a TU, grouped by operons
with open('tu_exp.tsv', 'w') as f:
	for genes in sets:
		tus = set(np.hstack(transcription.getTranscriptionUnitWithRna(list(genes))))
		singles = np.sum([1 if len(tu) == 1 else 0 for tu in transcription.getRnaInTranscriptionUnit(list(tus))])
		f.write('{} genes in {} TUs, with {} single genes\n'.format(len(genes), len(tus), singles))
		for tu in transcription.getRnaInTranscriptionUnit(list(tus)):
			f.write(json.dumps(tu.tolist()))
			for gene in tu:
				f.write('\t{}'.format(exp[gene]))
			f.write('\n')

# save counts for each operon
with open('gene_to_tu_count.tsv', 'w') as f:
	f.write('Gene\tTU\tSingle Genes\n')
	for genes in sets:
		tus = set(np.hstack(transcription.getTranscriptionUnitWithRna(list(genes))))
		singles = np.sum([1 if len(tu) == 1 else 0 for tu in transcription.getRnaInTranscriptionUnit(list(tus))])
		f.write('{}\t{}\t{}\n'.format(len(genes), len(tus), singles))

# calculate and save synthesis rates for each gene within an operon
with open('syn_rates.tsv', 'w') as f:
	f.write('Gene\tSynth Rate (exp * deg rate)\n')
	for genes in sets:
		tus = set(np.hstack(transcription.getTranscriptionUnitWithRna(list(genes))))
		singles = np.sum([1 if len(tu) == 1 else 0 for tu in transcription.getRnaInTranscriptionUnit(list(tus))])
		f.write('{} genes in {} TUs, with {} single genes\n'.format(len(genes), len(tus), singles))
		for tu in transcription.getRnaInTranscriptionUnit(list(tus)):
			f.write(json.dumps(tu.tolist()))
			for gene in tu:
				f.write('\t{}'.format(exp[gene]*deg[gene]))
			f.write('\n')
