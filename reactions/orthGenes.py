import numpy as np
import csv
import os
import cPickle
import re

# load sim data
fileLoc = os.path.dirname(os.path.dirname(__file__))
# sim_data = cPickle.load(open(os.path.join(fileLoc, 'scripts', 'sim_data.cp'), 'rb'))

# get genes from orth paper
orthGenes = []
reactionFile = csv.reader(open(os.path.join(fileLoc, 'reactions', 'orthReactions062317.tsv'), 'r'), delimiter='\t')
reactionFile.next()

# strip genes from the protein-reaction association column of the orth reactions
for line in reactionFile:
	genes = line[5]
	genes = re.sub('[\(\)]', '', genes)
	genes = re.sub(' and ', ' ', genes)
	genes = re.sub(' or ', ' ', genes)
	for gene in re.findall('[A-z]+', genes):
		if len(gene) > 4:
			if gene == 'SPONTANEOUS':
				continue

			# some have ec tacked to the end
			if gene[-2:] == 'ec':
				gene = gene[:-2]
			# multiple genes in same family are joined together eg TrpABC
			if len(gene) >= 4:
				for c in gene[3:]:
					orthGenes += [gene[:3] + c]
			else:
				orthGenes += [gene]
		elif len(gene) > 1:
			orthGenes += [gene]

orthGenes = np.unique(orthGenes)
orthGenes = [x.lower() for x in orthGenes]

# read information from formatted anotation file
annotatedGenes = []
annotationsFile = csv.reader(open(os.path.join(fileLoc, 'reactions', 'annotationsEdited062317.tsv'), 'r'), delimiter='\t')
annotationsFile.next()
annotationsFile.next()

for line in annotationsFile:
	annotatedGenes += [line[3].lower()]

import ipdb; ipdb.set_trace()