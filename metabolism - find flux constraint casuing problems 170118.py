#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- option to call a reduced form of metabolism (assume optimal)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

from itertools import izip

import numpy as np
from scipy.sparse import csr_matrix

import wholecell.processes.process
from wholecell.utils import units

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

from wholecell.utils.modular_fba import FluxBalanceAnalysis
from wholecell.utils.enzymeKinetics import EnzymeKinetics

from wholecell.utils.fitting import massesAndCountsToAddForHomeostaticTargets

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

SECRETION_PENALTY_COEFF = 1e-5

NONZERO_ENZYMES = False

USE_KINETICS = True
KINETICS_BURN_IN_PERIOD = 1

FBA_ITERATION_LIMIT = 100000
FBA_SOLVE_ITERATIONS = 5

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Metabolism, self).initialize(sim, sim_data)

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		self.exchangeConstraints = sim_data.process.metabolism.exchangeConstraints

		self.doublingTime = sim_data.doubling_time
		self.nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel

		self.getBiomassAsConcentrations = sim_data.mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel][0][1]
			)

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[sim_data.condition])
		concDict.update(self.concModificationsBasedOnCondition)

		self.objective = dict(
			(key, concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in concDict
			)

		self.getMass = sim_data.getter.getMass
		self.massReconstruction = sim_data.mass
		self.avgCellToInitialCellConvFactor = sim_data.mass.avgCellToInitialCellConvFactor

		self.ngam = sim_data.constants.nonGrowthAssociatedMaintenance

		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		self.energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		self.reactionStoich = sim_data.process.metabolism.reactionStoich
		self.externalExchangeMolecules = sim_data.nutrientData["secretionExchangeMolecules"]

		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel]:
			self.externalExchangeMolecules += sim_data.nutrientData["importExchangeMolecules"][nutrientsLabel]

			# Sorry
			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					nutrientsLabel, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)

		self.maintenanceReaction = sim_data.process.metabolism.maintenanceReaction
		self.externalExchangeMolecules = sorted(self.externalExchangeMolecules)
		self.extMoleculeMasses = self.getMass(self.externalExchangeMolecules)

		self.moleculeMasses = dict(zip(
			self.externalExchangeMolecules,
			self.getMass(self.externalExchangeMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)
			))

		# Data structures to compute reaction bounds based on enzyme presence/absence
		self.catalystsList = sim_data.process.metabolism.catalystsList
		self.reactionsWithCatalystsList = sim_data.process.metabolism.reactionCatalystsList
		self.reactionCatalystsDict = sim_data.process.metabolism.reactionCatalysts

		catalysisMatrixI = sim_data.process.metabolism.catalysisMatrixI
		catalysisMatrixJ = sim_data.process.metabolism.catalysisMatrixJ
		catalysisMatrixV = sim_data.process.metabolism.catalysisMatrixV

		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
		self.catalysisMatrix = csr_matrix((catalysisMatrixV, (catalysisMatrixI, catalysisMatrixJ)), shape = shape)

		self.catalyzedReactionBoundsPrev = np.inf * np.ones(len(self.reactionsWithCatalystsList))

		# Data structures to compute reaction targets based on kinetic parameters
		from reconstruction.ecoli.dataclasses.process.metabolism_constraints import constraints as kineticsConstraints
		self.kineticsConstraints = kineticsConstraints

		self.kineticsConstrainedReactions = sim_data.process.metabolism.constrainedReactionList
		# rxns = ['ALARACECAT-RXN (reverse)','ALARACECAT-RXN', 'ASNSYNA-RXN', 'ASPCARBTRANS-RXN', 'DIHYDROOROT-RXN (reverse)', 'GLUTAMINESYN-RXN', 'IMP-DEHYDROG-RXN', 'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.', 'OROTPDECARB-RXN', 'PHOSACETYLTRANS-RXN (reverse)', 'PSERTRANSAM-RXN', 'RIBOSYLHOMOCYSTEINASE-RXN', 'RXN-11832', 'RXN-14325', 'RXN-15943', 'RXN0-1461', 'RXN0-5114[CCO-CYTOSOL]-3-P-SERINE/WATER//SER/Pi.38.', 'RXN0-7008-PRO/UBIQUINONE-8//L-DELTA1-PYRROLINE_5-CARBOXYLATE/CPD-9956/PROTON.67.', 'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.', 'URPHOS-RXN (reverse)']
		# rxns = ['ASNSYNA-RXN', 'ASPCARBTRANS-RXN', 'DIHYDROOROT-RXN (reverse)', 'GLUTAMINESYN-RXN', 'IMP-DEHYDROG-RXN', 'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.', 'OROTPDECARB-RXN', 'PHOSACETYLTRANS-RXN (reverse)', 'PSERTRANSAM-RXN', 'RXN-11832', 'RXN0-5114[CCO-CYTOSOL]-3-P-SERINE/WATER//SER/Pi.38.', 'RXN0-7008-PRO/UBIQUINONE-8//L-DELTA1-PYRROLINE_5-CARBOXYLATE/CPD-9956/PROTON.67.', 'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.']
		
		self.excluded = []
		# self.included = np.array(self.kineticsConstrainedReactions, dtype = bool)
		# for rxn in rxns[::-1]:
		# 	self.excluded.append(self.kineticsConstrainedReactions.index(rxn))
		# for idx in self.excluded:
		# 	self.included[idx] = False
		# 	self.kineticsConstrainedReactions.pop(idx)
		self.kineticsEnzymesList = sim_data.process.metabolism.enzymeIdList
		self.kineticsSubstratesList = sim_data.process.metabolism.kineticsSubstratesList

		constraintToReactionMatrixI = sim_data.process.metabolism.constraintToReactionMatrixI
		constraintToReactionMatrixJ = sim_data.process.metabolism.constraintToReactionMatrixJ
		constraintToReactionMatrixV = sim_data.process.metabolism.constraintToReactionMatrixV
		shape = (constraintToReactionMatrixI.max() + 1, constraintToReactionMatrixJ.max() + 1)
		self.constraintToReactionMatrix = np.zeros(shape, np.float64)
		self.constraintToReactionMatrix[constraintToReactionMatrixI, constraintToReactionMatrixJ] = constraintToReactionMatrixV
		self.constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly

		self.metabolismKineticObjectiveWeight = sim_data.constants.metabolismKineticObjectiveWeight

		# Set up FBA solver
		self.fbaObjectOptions = {
			"reactionStoich" : self.reactionStoich,
			"externalExchangedMolecules" : self.externalExchangeMolecules,
			"objective" : self.objective,
			"objectiveType" : "homeostatic_kinetics_mixed",
			"objectiveParameters" : {
					"kineticObjectiveWeight":1e-7,#self.metabolismKineticObjectiveWeight,
					"reactionRateTargets":{reaction: 1 for reaction in self.kineticsConstrainedReactions}, #This target is arbitrary, it gets reset each timestep during evolveState
					"oneSidedReactionTargets":[], # TODO: deal with this dynamically (each time step)
					},
			"moleculeMasses" : self.moleculeMasses,
			"secretionPenaltyCoeff" : SECRETION_PENALTY_COEFF, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver" : "glpk",
			"maintenanceCostGAM" : self.energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction" : self.maintenanceReaction,
		}

		if USE_KINETICS == False:
			self.fbaObjectOptions["objectiveType"] = "homeostatic"

		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
		self.fba._solver._model.set_iteration_limit(FBA_ITERATION_LIMIT)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])

		# Disable all rates during burn-in
		if USE_KINETICS and KINETICS_BURN_IN_PERIOD > 0:
			self.fba.disableKineticTargets()
			self.burnInComplete = False
		else:
			self.burnInComplete = True

		self.currentNgam = 1 * (COUNTS_UNITS / VOLUME_UNITS)
		self.currentPolypeptideElongationEnergy = 1 * (COUNTS_UNITS / VOLUME_UNITS)

		# Set constraints
		## External molecules
		self.externalMoleculeIDs = self.fba.getExternalMoleculeIDs()

		# Views
		self.metaboliteNames = self.fba.getOutputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.catalystsList)
		self.kineticsEnzymes = self.bulkMoleculesView(self.kineticsEnzymesList)
		self.kineticsSubstrates = self.bulkMoleculesView(self.kineticsSubstratesList)

		outputMoleculeIDs = self.fba.getOutputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.getInternalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.catalysts.requestAll()
		self.kineticsEnzymes.requestAll()
		self.kineticsSubstrates.requestAll()

	def evolveState(self):
		# Solve for metabolic fluxes
		metaboliteCountsInit = self.metabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

		cellVolume = cellMass / self.cellDensity

		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		self.newPolypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			self.newPolypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(self.nutrientToDoublingTime.get(self._sim.processes["PolypeptideElongation"].currentNutrients, self.nutrientToDoublingTime["minimal"]))

		# Set external molecule levels
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)

		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			self.nutrientsTimeSeriesLabel,
			self.time(),
			self.concModificationsBasedOnCondition,
			)

		if newObjective != None and newObjective != self.objective:
			# Build new fba instance with new objective
			self.fbaObjectOptions["objective"] = newObjective
			self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
			self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])

		# After completing the burn-in, enable kinetic rates
		if (USE_KINETICS) and (not self.burnInComplete) and (self._sim.time() > KINETICS_BURN_IN_PERIOD):
			self.burnInComplete = True
			self.fba.enableKineticTargets()

		# Set external molecule levels
		self.fba.setExternalMoleculeLevels(externalMoleculeLevels)

		self.newNgam = self.ngam * coefficient

		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		ADJUSTMENT_RATIO = .01
		ngam_diff = np.abs(self.currentNgam.asNumber() - self.newNgam.asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		if ngam_diff > ADJUSTMENT_RATIO:
			self.currentNgam = self.newNgam
			self.fba.setMaxReactionFlux(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.fba.setMinReactionFlux(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))

		poly_diff = np.abs((self.currentPolypeptideElongationEnergy.asNumber() - self.newPolypeptideElongationEnergy.asNumber())) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		if poly_diff > ADJUSTMENT_RATIO:
			self.currentPolypeptideElongationEnergy = self.newPolypeptideElongationEnergy
			self.fba.setMaxReactionFlux(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.fba.setMinReactionFlux(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		#  Find metabolite concentrations from metabolite counts
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit[self.internalExchangeIdxs]

		# Make a dictionary of metabolite names to metabolite concentrations
		metaboliteConcentrationsDict = dict(zip(self.metaboliteNames, metaboliteConcentrations))

		self.fba.setInternalMoleculeLevels(
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS)
			)

		catalystsCountsInit = self.catalysts.counts()

		kineticsEnzymesCountsInit = self.kineticsEnzymes.counts()
		kineticsEnzymesConcentrations = countsToMolar * kineticsEnzymesCountsInit

		kineticsSubstratesCountsInit = self.kineticsSubstrates.counts()
		kineticsSubstratesConcentrations = countsToMolar * kineticsSubstratesCountsInit

		if NONZERO_ENZYMES:
			# Add one of every enzyme to ensure none are zero
			catalystsCountsInit += 1
			kineticsEnzymesConcentrations = countsToMolar * (kineticsEnzymesCountsInit + 1)

		# This should kill the cell
		# catalystsCountsInit[136] *= 0
		# catalystsCountsInit[940] *= 0

		if USE_KINETICS and self.burnInComplete:

			# Set hard upper bounds constraints based on enzyme presence (infinite upper bound) or absence (upper bound of zero)
			catalyzedReactionBounds = np.inf * np.ones(len(self.reactionsWithCatalystsList))
			rxnPresence = self.catalysisMatrix.dot(catalystsCountsInit)
			catalyzedReactionBounds[rxnPresence == 0] = 0

			updateIdxs = np.where(catalyzedReactionBounds != self.catalyzedReactionBoundsPrev)[0]
			updateRxns = [self.reactionsWithCatalystsList[idx] for idx in updateIdxs]
			updateVals = catalyzedReactionBounds[updateIdxs]

			self.fba.setMaxReactionFluxes(updateRxns, updateVals, raiseForReversible = False)

			self.catalyzedReactionBoundsPrev = catalyzedReactionBounds

			# Set target fluxes for reactions based on their most relaxed constraint
			constraintValues = self.kineticsConstraints(
				kineticsEnzymesConcentrations.asNumber(units.umol / units.L),
				kineticsSubstratesConcentrations.asNumber(units.umol / units.L),
				)
			reactionTargets = (units.umol / units.L / units.s) * np.max(self.constraintToReactionMatrix * constraintValues, axis = 1)
			# import ipdb; ipdb.set_trace()
			# record which constraint was used, add constraintToReactionMatrix to ensure the index is one of the constraints if multiplication is 0
			reactionConstraint = np.argmax(self.constraintToReactionMatrix * constraintValues + self.constraintToReactionMatrix, axis = 1)

			targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(COUNTS_UNITS / VOLUME_UNITS)
			targets = np.delete(targets, self.excluded)
			self.fba.setKineticTarget(
				self.kineticsConstrainedReactions,
				targets,
				raiseForReversible = False,
				)
			# self.fba.disableKineticTargets()
			# import ipdb; ipdb.set_trace()
			# self.setCentralCarbonMetabolismFluxes(coefficient)
			self.fba.disableKineticTargets("SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.")
			self.fba.disableKineticTargets("RXN0-7008-PRO/UBIQUINONE-8//L-DELTA1-PYRROLINE_5-CARBOXYLATE/CPD-9956/PROTON.67.")
			self.fba.disableKineticTargets("R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.")
			self.fba.disableKineticTargets("FUMHYDR-RXN")
			# self.fba.disableKineticTargets("GLUTAMINESYN-RXN")
			# self.fba.disableKineticTargets("GLUTDEHYD-RXN")
			# self.fba.disableKineticTargets("ASNSYNA-RXN")
			# self.fba.disableKineticTargets("DIHYDROOROT-RXN (reverse)")
			# self.fba.disableKineticTargets("PRIBFAICARPISOM-RXN")
			# self.fba.disableKineticTargets("PSERTRANSAM-RXN")
			# self.fba.disableKineticTargets("ALARACECAT-RXN (reverse)")
			# import ipdb; ipdb.set_trace()
			# if self._sim.time() > 10:
			# b = True
			# for rxn in self.kineticsConstrainedReactions:
			# 	# if rxn == "PRIBFAICARPISOM-RXN" or rxn == "PSERTRANSAM-RXN" or rxn == "ASNSYNA-RXN" or rxn == "DIHYDROOROT-RXN (reverse)" or rxn == "FUMHYDR-RXN" or rxn == "GLUTAMINESYN-RXN" or rxn == "GLUTDEHYD-RXN" or rxn == "RXN0-7008-PRO/UBIQUINONE-8//L-DELTA1-PYRROLINE_5-CARBOXYLATE/CPD-9956/PROTON.67." or rxn == "SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.":
			# 	# 	continue
			# 	if rxn == "SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31." or rxn == "RXN0-7008-PRO/UBIQUINONE-8//L-DELTA1-PYRROLINE_5-CARBOXYLATE/CPD-9956/PROTON.67.":
			# 		continue
			# 	if rxn == "GLUTRACE-RXN (reverse)" or rxn == "GLUTRNAREDUCT-RXN-GLUTAMATE-1-SEMIALDEHYDE/GLT-tRNAs/NADP//Charged-GLT-tRNAs/NADPH.65.":
			# 		b = False
			# 		continue
			# 	if b:
			# 		continue
			# 	# if rxn == "FUMHYDR-RXN" or rxn == "ALARACECAT-RXN (reverse)":
			# 	# 	continue
			# 	self.fba.disableKineticTargets(rxn)
			# 	glc = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.externalExchangeFluxes()[31] / coefficient).asNumber(units.mmol / units.g / units.h)
			# 	if glc < -6.54:
			# 		print "** %s, %f" % (rxn, glc)
			# 	else:
			# 		print rxn, glc
			# 	self.fba.enableKineticTargets(rxn)
			# import ipdb; ipdb.set_trace()



		# if self._sim.time() > 0.8:
		# 	for reaction in self.kineticsConstrainedReactions:
		# 		self.fbaObjectOptions = {
		# 			"reactionStoich" : self.reactionStoich,
		# 			"externalExchangedMolecules" : self.externalExchangeMolecules,
		# 			"objective" : self.objective,
		# 			"objectiveType" : "homeostatic_kinetics_mixed",
		# 			"objectiveParameters" : {
		# 					"kineticObjectiveWeight":1e-7,#self.metabolismKineticObjectiveWeight,
		# 					"reactionRateTargets":{reaction: 1e-5}, #This target is arbitrary, it gets reset each timestep during evolveState
		# 					"oneSidedReactionTargets":[], # TODO: deal with this dynamically (each time step)
		# 					},
		# 			"moleculeMasses" : self.moleculeMasses,
		# 			"secretionPenaltyCoeff" : SECRETION_PENALTY_COEFF, # The "inconvenient constant"--limit secretion (e.g., of CO2)
		# 			"solver" : "glpk",
		# 			"maintenanceCostGAM" : self.energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
		# 			"maintenanceReaction" : self.maintenanceReaction,
		# 		}
		# 		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
		# 		self.fba._solver._model.set_iteration_limit(FBA_ITERATION_LIMIT)
		# 		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		# 		self.newNgam = self.ngam * coefficient

		# 		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		# 		ADJUSTMENT_RATIO = .01
		# 		ngam_diff = np.abs(self.currentNgam.asNumber() - self.newNgam.asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		# 		if ngam_diff > ADJUSTMENT_RATIO:
		# 			self.currentNgam = self.newNgam
		# 			self.fba.maxReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))
		# 			self.fba.minReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# 		poly_diff = np.abs((self.currentPolypeptideElongationEnergy.asNumber() - self.newPolypeptideElongationEnergy.asNumber())) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		# 		if poly_diff > ADJUSTMENT_RATIO:
		# 			self.currentPolypeptideElongationEnergy = self.newPolypeptideElongationEnergy
		# 			self.fba.maxReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		# 			self.fba.minReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# 		#  Find metabolite concentrations from metabolite counts
		# 		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit[self.internalExchangeIdxs]

		# 		# Make a dictionary of metabolite names to metabolite concentrations
		# 		metaboliteConcentrationsDict = dict(zip(self.metaboliteNames, metaboliteConcentrations))

		# 		self.fba.internalMoleculeLevelsIs(
		# 			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS)
		# 			)

		# 		catalystsCountsInit = self.catalysts.counts()

		# 		kineticsEnzymesCountsInit = self.kineticsEnzymes.counts()
		# 		kineticsEnzymesConcentrations = countsToMolar * kineticsEnzymesCountsInit

		# 		kineticsSubstratesCountsInit = self.kineticsSubstrates.counts()
		# 		kineticsSubstratesConcentrations = countsToMolar * kineticsSubstratesCountsInit

		# 		if NONZERO_ENZYMES:
		# 			# Add one of every enzyme to ensure none are zero
		# 			catalystsCountsInit += 1
		# 			kineticsEnzymesConcentrations = countsToMolar * (kineticsEnzymesCountsInit + 1)

		# 		# This should kill the cell
		# 		# catalystsCountsInit[136] *= 0
		# 		# catalystsCountsInit[940] *= 0

		# 		if USE_KINETICS and self.burnInComplete:

		# 			# Set hard upper bounds constraints based on enzyme presence (infinite upper bound) or absence (upper bound of zero)
		# 			catalyzedReactionBounds = np.inf * np.ones(len(self.reactionsWithCatalystsList))
		# 			rxnPresence = self.catalysisMatrix.dot(catalystsCountsInit)
		# 			catalyzedReactionBounds[rxnPresence == 0] = 0

		# 			updateIdxs = np.where(catalyzedReactionBounds != self.catalyzedReactionBoundsPrev)[0]
		# 			updateRxns = [self.reactionsWithCatalystsList[idx] for idx in updateIdxs]
		# 			updateVals = catalyzedReactionBounds[updateIdxs]

		# 			self.fba.setMaxReactionFluxes(updateRxns, updateVals, raiseForReversible = False)

		# 			self.catalyzedReactionBoundsPrev = catalyzedReactionBounds

		# 			# Set target fluxes for reactions based on their most relaxed constraint
		# 			constraintValues = self.kineticsConstraints(
		# 				kineticsEnzymesConcentrations.asNumber(units.umol / units.L),
		# 				kineticsSubstratesConcentrations.asNumber(units.umol / units.L),
		# 				)
		# 			reactionTargets = (units.umol / units.L / units.s) * np.max(self.constraintToReactionMatrix * constraintValues, axis = 1)

		# 			# record which constraint was used, add constraintToReactionMatrix to ensure the index is one of the constraints if multiplication is 0
		# 			reactionConstraint = np.argmax(self.constraintToReactionMatrix * constraintValues + self.constraintToReactionMatrix, axis = 1)

		# 			targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(COUNTS_UNITS / VOLUME_UNITS)
		# 			import ipdb; ipdb.set_trace()
		# 			self.fba.setKineticTarget(
		# 				self.kineticsConstrainedReactions,
		# 				(TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(COUNTS_UNITS / VOLUME_UNITS)[self.included],
		# 				raiseForReversible = False,
		# 				)
		# 			# import ipdb; ipdb.set_trace()
		# 			# self.setCentralCarbonMetabolismFluxes(coefficient)
		# 			# self.fba.disableKineticTargets("SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.")
		# 			# self.fba.disableKineticTargets("R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.")
		# 			self.fba.disableKineticTargets()
		# 			# import ipdb; ipdb.set_trace()
		# 			# if self._sim.time() > 10:
		# 			# for rxn in self.kineticsConstrainedReactions:
		# 			# 	# if rxn == "SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.":
		# 			# 	# 	continue
		# 			# 	self.fba.enableKineticTargets(rxn)
		# 			# 	if self.fba.externalExchangeFluxes()[31] < -0.0045:
		# 			# 		print "** %s, %f" % (rxn, self.fba.externalExchangeFluxes()[31])
		# 			# 	else:
		# 			# 		print rxn, self.fba.externalExchangeFluxes()[31]
		# 			# 	self.fba.disableKineticTargets(rxn)
		# 			# import ipdb; ipdb.set_trace()
		# 		glc = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.externalExchangeFluxes()[31] / coefficient).asNumber(units.mmol / units.g / units.h)
		# 		if glc < -5.161:
		# 			print reaction, glc
		# 	import ipdb; ipdb.set_trace()

		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.getOutputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)
		if USE_KINETICS and self.burnInComplete:
			relError = np.abs((self.fba.getReactionFluxes(self.kineticsConstrainedReactions) - targets) / (targets + 1e-15))


		exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.getExternalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)
		print exFluxes[31]

		# TODO: report as reactions (#) per second & store volume elsewhere

		self.writeToListener("FBAResults", "deltaMetabolites",
			deltaMetabolites)

		self.writeToListener("FBAResults", "reactionFluxes",
							 self.fba.getReactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes",
			exFluxes)

		self.writeToListener("FBAResults", "objectiveValue",  # TODO
							 self.fba.getObjectiveValue()) # / len(deltaMetabolites)) # divide to normalize by number of metabolites

		self.writeToListener("FBAResults", "rowDualValues",
							 self.fba.getShadowPrices(self.metaboliteNames))

		self.writeToListener("FBAResults", "columnDualValues",
							 self.fba.getReducedCosts(self.fba.getReactionIDs()))

		# self.writeToListener("FBAResults", "kineticObjectiveValues",
		# 	self.fba.kineticObjectiveValues())

		# self.writeToListener("FBAResults", "homeostaticObjectiveValues",
		# 	self.fba.homeostaticObjectiveValues())

		# self.writeToListener("FBAResults", "homeostaticObjectiveWeight",
		# 	self.fba.homeostaticObjectiveWeight())

		# self.writeToListener("EnzymeKinetics", "baseRates",
		# 	self.baseRates.asNumber(FLUX_UNITS))

		# self.writeToListener("EnzymeKinetics", "reactionKineticPredictions",
		# 	self.allRateEstimates.asNumber(FLUX_UNITS))

		# self.writeToListener("EnzymeKinetics", "overconstraintMultiples",
		# 	self.overconstraintMultiples)

		# self.writeToListener("EnzymeKinetics", "kineticTargetFluxes",
		# 	self.fba.kineticTargetFluxes())

		# self.writeToListener("EnzymeKinetics", "kineticTargetErrors",
		# 	self.fba.kineticTargetFluxErrors())

		# self.writeToListener("EnzymeKinetics", "kineticTargetRelativeDifferences",
		# 	self.fba.kineticTargetFluxRelativeDifferences())

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit",
			metaboliteCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal",
			metaboliteCountsFinal)

		self.writeToListener("EnzymeKinetics", "enzymeCountsInit",
			kineticsEnzymesCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations",
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		self.writeToListener("EnzymeKinetics", "countsToMolar",
			countsToMolar.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		self.writeToListener("EnzymeKinetics", "actualFluxes",
							 self.fba.getReactionFluxes(self.kineticsConstrainedReactions) / self.timeStepSec())

		if USE_KINETICS and self.burnInComplete:
			self.writeToListener("EnzymeKinetics", "targetFluxes",
				targets / self.timeStepSec())

			self.writeToListener("EnzymeKinetics", "reactionConstraint",
				reactionConstraint)

	def setCentralCarbonMetabolismFluxes(self, coefficient):
		toyaFluxes = {"TRANS-RXN-157": 10.17,
			"TRANS-RXN-157 (reverse)": 0,
			"RXN0-7077-GLC/PROTON//GLC/PROTON.23.": 10.17,
			"RXN0-7077-GLC/PROTON//GLC/PROTON.23. (reverse)": 0,
			"PGLUCISOM-RXN": 0,
			"PGLUCISOM-RXN (reverse)": 6.04,
			"6PFRUCTPHOS-RXN": 8.06,
			"6PFRUCTPHOS-RXN (reverse)": 0,
			"F16ALDOLASE-RXN": 0,
			"F16ALDOLASE-RXN (reverse)": 8.06,
			"TRIOSEPISOMERIZATION-RXN": 8.06,
			"TRIOSEPISOMERIZATION-RXN (reverse)": 0,
			"2PGADEHYDRAT-RXN": 0,
			"2PGADEHYDRAT-RXN (reverse)": 16.89,
			"RXN-15513": 0,
			"RXN-15513 (reverse)": 15.85,
			"3PGAREARR-RXN": 15.85,
			"3PGAREARR-RXN (reverse)": 0,
			"PEPDEPHOS-RXN": 0,
			"PEPDEPHOS-RXN (reverse)": 1.97,
			"PYRUVDEH-RXN": 11.12,
			"PYRUVDEH-RXN (reverse)": 0,
			"GLU6PDEHYDROG-RXN": 3.98,
			"GLU6PDEHYDROG-RXN (reverse)": 0,
			"RXN-9952": 3.98,
			"RXN-9952 (reverse)": 0,
			"RIBULP3EPIM-RXN": 2.09,
			"RIBULP3EPIM-RXN (reverse)": 0,
			"RIB5PISOM-RXN": 0,
			"RIB5PISOM-RXN (reverse)": 1.88,
			"1TRANSKETO-RXN": 0,
			"1TRANSKETO-RXN (reverse)": 1.20,
			"TRANSALDOL-RXN": 0,
			"TRANSALDOL-RXN (reverse)": 1.20,
			"2TRANSKETO-RXN": 0.89,
			"2TRANSKETO-RXN (reverse)": 0,
			"CITSYN-RXN": 5.46,
			"CITSYN-RXN (reverse)": 0,
			"ISOCITDEH-RXN": 5.46,
			"ISOCITDEH-RXN (reverse)": 0,
			"2OXOGLUTARATEDEH-RXN": 4.62,
			"2OXOGLUTARATEDEH-RXN (reverse)": 0,
			"SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.": 4.62,
			"SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31. (reverse)": 0,
			"FUMHYDR-RXN": 4.62,
			"FUMHYDR-RXN (reverse)": 0,
			"MALATE-DEH-RXN": 3.50,
			"MALATE-DEH-RXN (reverse)": 0,
			"PEPCARBOX-RXN": 3.34,
			"PEPCARBOX-RXN (reverse)": 0,
			"MALIC-NADP-RXN": 1.12,
			"MALIC-NADP-RXN (reverse)": 0,
			}

		fbaReactions = self.fba.getReactionIDs()
		for reaction in toyaFluxes:
			if reaction in fbaReactions:
				flux = ((units.mmol / units.g / units.h) * coefficient * toyaFluxes[reaction]).asNumber(COUNTS_UNITS / VOLUME_UNITS)
				self.fba.setMinReactionFlux(reaction, 0.3 * flux, raiseForReversible = False)
				self.fba.setMaxReactionFlux(reaction, 1.7 * flux, raiseForReversible = False)
