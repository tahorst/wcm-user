from models.ecoli.sim.simulation import EcoliDaughterSimulation
import os

j = 0	# seed in N_INIT_SIMS
k = 2	# gen
l = 0	# cell, 0 if single daugther

path = "/home/thorst/wcEcoli/out/SET_D/nutrientTimeSeries_000002/"

options = {}		
options["simDataLocation"] = os.path.join(path, "kb/simData_Modified.cPickle")
options["inheritedStatePath"] = os.path.join(path, "000000/generation_000001/000000/simOut/Daughter1")
options["logToDisk"] = False
options["overwriteExistingFiles"] = False
options["seed"] = (j + 1) * ((2**k - 1) + l)
options["massDistribution"] = 1
options["growthRateNoise"] = 1
options["dPeriodDivision"] = 1
options["translationSupply"] = 1

sim = EcoliDaughterSimulation(**options)

# sim = EcoliDaughterSimulation(simDataLocation = os.path.join(path, "kb/simData_Modified.cPickle"),
# 	inheritedStatePath = os.path.join(path, "000044/generation_000000/000000/simOut/Daughter1"),
# 	overwriteExistingFiles = False,
# 	seed = (j + 1) * ((2**k - 1) + l),
# 	mass_distribution = 1,
# 	growth_rate_noise = 1,
# 	d_period_division = 1,
# 	translation_supply = 1)
sim.run()