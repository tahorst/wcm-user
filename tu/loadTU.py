'''
Produces and saves a dict for tu as json object in OUTPUT based on annotations from Mialy
Uses sequencing data set to identify which genes form TUs
Loads information from MATRIX, ROWS and COLS files
Stores matrices by indices for more efficient storage because they are sparse
A is matrix mapping genes to transcription units (genes by TU)
dict contains keys:
	mati, matj - indices for A matrix (values of 1)
	pinvi, pinvj, pinvv - indices and values for pinv of A
	dims - dimensions of A matrix
	rows, cols - names of genes and list of genes within a TU corresponding to A
'''

import json
import cPickle
import numpy as np
import os

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

# file references
DIR = '/home/thorst/wcEcoli/user/tu'
MATRIX = 'geneTUMatrix.pkl'
ROWS = 'geneTUMatrix_Rows.pkl'
COLS = 'geneTUMatrix_Columns.pkl'
OUTPUT = 'tu.json'


raw_data = KnowledgeBaseEcoli()

# read cPickle data from Mialy
with open(os.path.join(DIR, MATRIX)) as f:
	mat = np.array(cPickle.load(f))
with open(os.path.join(DIR, ROWS)) as f:
	rows = cPickle.load(f)
with open(os.path.join(DIR, COLS)) as f:
	cols = cPickle.load(f)

# map ecocyc ID to wcm ID for row and col annotation
ecToWcGeneMap =  {gene['id']: gene['rnaId'] for gene in raw_data.genes}
wcmRows = [ecToWcGeneMap[gene] for gene in rows]
wcmCols = [[ecToWcGeneMap.get(gene, '') for gene in col] for col in cols]

# analyze matrix
print '%i genes without an transcription unit' % (np.sum(np.sum(mat, axis=1) == 0))
print '%i transcription units without a gene' % (np.sum(np.sum(mat, axis=0) == 0))

# represent matrix more compactly
mati, matj = np.where(mat)

# create new entries mapping genes that do not appear in data to their own operon
for geneIdx in np.where(np.sum(mat, axis=1) == 0)[0]:
	mati = np.hstack((mati, geneIdx))
	matj = np.hstack((matj, len(wcmCols)))
	wcmCols.append([wcmRows[geneIdx]])

# sort rows by name to match order in sim_data for synth prob
sortedIndices = {j: i for i, j in enumerate(np.argsort(wcmRows))}
matiSorted = np.array([sortedIndices[idx] for idx in mati])

# reconstruct matrix with new columns
dims = (len(wcmRows), len(wcmCols))
expandedMat = np.zeros(dims)
expandedMat[matiSorted, matj] = 1

# create pseudoinverse for least norm analysis
## TODO - is removal of low values justified?
## TODO - explicit least norm formulation instead of pinv - use SVD
pinv = np.linalg.pinv(expandedMat)
pinvi, pinvj = np.where(np.abs(pinv) > 1e-12)
pinvv = pinv[pinvi, pinvj]

# create dict with all information
tu = {}
tu['mati'] = matiSorted.tolist()
tu['matj'] = matj.tolist()
tu['rows'] = sorted(wcmRows)
tu['cols'] = wcmCols
tu['dims'] = dims
tu['pinvi'] = pinvi.tolist()
tu['pinvj'] = pinvj.tolist()
tu['pinvv'] = pinvv.tolist()

import ipdb; ipdb.set_trace()

# save json object
with open(os.path.join(DIR, OUTPUT), 'w') as f:
	json.dump(tu, f)
