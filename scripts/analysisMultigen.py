# runs all multigen analysis plots for a given sim
# if no folder is given it runs for most recent in wcEcoli_out other picks first arg passed
# only need parent directory in wcEcoli_out as argument

import os
import sys
import cPickle

from wholecell.fireworks.firetasks.analysisMultiGen import AnalysisMultiGenTask

BASE_DIR = "/scratch/users/thorst/wcEcoli_out/"
SEED = "000000"

if len(sys.argv) < 2:
	for folder in sorted(os.listdir(BASE_DIR), reverse = True):
		if folder[0].isdigit():
			path = os.path.join(BASE_DIR, folder)
			break
else:
	arg = sys.argv[1]
	if arg.startswith("out/"):
		arg = arg[4:]
	path = os.path.join(BASE_DIR, arg)

for folder in os.listdir(path):
	if "_" in folder:
		variant = folder

# variant = "condition_000000"
# variant = "geneKnockdown_000030"

inputSeedDirectory = os.path.join(path, variant, SEED)
outputDir = os.path.join(path, variant, SEED, "plotOut")
simData = os.path.join(path, variant, "kb/simData_Modified.cPickle")
validationData = os.path.join(path, "kb/validationData.cPickle")
with open(os.path.join(path, 'metadata', 'metadata.cPickle')) as f:
	metadata = cPickle.load(f)
metadata['analysis_type'] = 'multigen'
metadata['variant_function'] = variant
metadata['variant_index'] = None
metadata['seed'] = SEED

print('Multigen analysis from {}'.format(inputSeedDirectory))
task = AnalysisMultiGenTask(input_seed_directory=inputSeedDirectory, input_sim_data=simData, input_validation_data=validationData, output_plots_directory=outputDir, metadata=metadata)
task.run_task(None)
