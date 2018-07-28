'''
Derivation of growth control loosely based on Bosdriesz et al
Addition of protein, mRNA, and more comprehensive growth rate

TODO - realistic values for AA (fraction in protein and synthesis rates)
TODO - dependent on SpoT?

Not working - can't get protein concentrations to be handled properly
Upshift appears to be close to working
May need to double check rates to make sure using concentrations
'''

from __future__ import division

from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

# indices for concentrations
n_aa = 20
n_species = 6  # species of RNA and protein (RelA, SpoT, RNAP, ribosome, metabolism, synthetases)
RELA = 0
SPOT = 1
RNAP = 2
RIBOSOME = 3
METABOLISM = 4
SYNTHETASE = 5
aa_indices = range(n_aa)
trna_indices = range(n_aa, 2*n_aa)
trna_charged_indices = range(2*n_aa, 3*n_aa)
start = trna_charged_indices[-1] + 1
mrna_indices = range(start, start + n_species)
start = mrna_indices[-1] + 1
protein_indices = range(start, start + n_species)
ppgpp_index = protein_indices[-1] + 1
rnap_bound_index = ppgpp_index + 1
biomass_index = rnap_bound_index + 1

n_molecules = biomass_index + 1

# parameters
## amino acids
k_aa_syn = 0.2 * np.random.lognormal(np.log(1), 0.2, 20)  # AA production rate, from text under S4 equations and mathematica file (Bosdriesz)
KI_aa = 100  # uM - inhibition constant for synthesis (Bosdriesz)
f_aa = 1 / 20  # fraction of AA in proteins that are a given species

## charging
k_synthetase = 100  # 1/s - kcat for synthetase
KM_trna_uncharged = 1  # uM - Michaelis constant for synthetases for uncharged tRNA
KM_aa = 100  # uM - Michaelis constant for synthetases for AA

## ribosome
rho = 1 / 20
k_rib_elong = 20  # 1/s - protein elongation rate
KD_trna_charged = 1  # uM - charged tRNA-ribosome dissociation constant (Bosdriesz)
KD_trna_uncharged = 500  # uM - uncharged tRNA-ribosome dissociation constant (Bosdriesz)

## ppGpp
k_rela = 75  # 1/s - RelA ppGpp synthesis (Bosdriesz)
KD_rela = 0.26  # uM - uncharged tRNA conc for half maximal ppGpp production (Bosdriesz)
v_spot_syn = 0.001  # uM/s - SpoT ppGpp synthesis (Bosdriesz)
k_spot_deg = np.log(2) / 30  # 1/s - ppGpp half life (Bosdriesz)

## RNAP - TODO - get rates and KI
k_rnap_elong = 50  # nt/s - RNAP elongation rate
k_bind = 0.2
KM_rnap_free = 20  # uM - Michaelis constant for rrn-promoters (Bosdriesz)
k_trna = k_rnap_elong / 80  # 1/s - tRNA length is ~80 nts
KI_trna = 1  # uM - ppGpp inhibition of tRNA transcription
protein_lengths = np.array([
	744,
	702,
	3400,
	7000,
	300,
	550,
	])
rna_lengths = protein_lengths * 3
k_mrna = k_rnap_elong / rna_lengths
KI_mrna = np.array([
	10,
	10,
	5,
	1,
	10,
	5,
	])

## additional parameters


def dcdt(c, t):
	# amino acid synthesis shift
	nutrients = 1
	if t > 400:
		nutrients = 0.1
	if t > 400:
		nutrients = 10

	dc = np.zeros_like(c)

	# concentrations
	aa = c[aa_indices]
	trna = c[trna_indices]
	trna_charged = c[trna_charged_indices]
	mrna = c[mrna_indices]
	protein = c[protein_indices]
	ppgpp = c[ppgpp_index]
	rnap_bound = c[rnap_bound_index]
	biomass = c[biomass_index]

	# derived concentrations
	## proteins
	rela = protein[RELA]
	spot = protein[SPOT]
	rnap = protein[RNAP]
	ribosome = protein[RIBOSOME]
	metabolism = protein[METABOLISM]
	synthetase = protein[SYNTHETASE]

	## rates/parameters
	k_trna_syn = k_trna / (1 + ppgpp / KI_trna)
	k_mrna_syn = k_mrna / (1 + ppgpp / KI_mrna)
	k_completion = k_trna_syn + np.sum(k_mrna_syn)  # rate of RNAP completion, needs to be < 1
	f_protein = mrna / np.sum(mrna)

	## dependent species
	trna_uncharged = trna - trna_charged
	ribosome_uncharged = ribosome * rho * trna_uncharged / KD_trna_uncharged / (1 + trna_charged / KD_trna_charged + trna_uncharged / KD_trna_uncharged)
	ribosome_uncharged_total = np.sum(ribosome_uncharged)
	rnap_free = rnap - rnap_bound

	# rates
	v_aa_syn = nutrients * k_aa_syn * metabolism / (1 + aa / KI_aa)
	v_trna_syn = k_trna_syn * rnap_bound
	v_charging = k_synthetase * synthetase * trna_uncharged / KM_trna_uncharged * aa / KM_aa / (1 + trna_uncharged / KM_trna_uncharged + aa / KM_aa + trna_uncharged / KM_trna_uncharged * aa / KM_aa)
	v_ribosome = k_rib_elong * ribosome / (1 + f_aa * np.sum(KD_trna_charged / trna_charged + KD_trna_charged / trna_charged * trna_uncharged / KD_trna_uncharged))
	v_transcription = k_mrna_syn * rnap_bound
	v_rela = k_rela * rela * ribosome_uncharged_total / KD_rela / (1 + ribosome_uncharged_total / KD_rela)
	# v_spot_syn = 1
	v_spot_deg = k_spot_deg * ppgpp
	v_rnap_binding = k_bind * rnap_free / (KM_rnap_free + rnap_free)
	v_rnap_completion = k_completion * rnap_bound

	mu = (np.sum(v_aa_syn/aa) + np.sum(v_trna_syn/trna) + np.sum(v_transcription/mrna)) / biomass

	# derivatives
	dc[aa_indices] = v_aa_syn - v_charging - mu * aa
	dc[trna_indices] = v_trna_syn - mu * trna
	dc[trna_charged_indices] = v_charging - f_aa * v_ribosome - mu * trna_charged
	dc[mrna_indices] = v_transcription - mu * mrna
	dc[protein_indices] = f_protein * v_ribosome / protein_lengths - mu * protein
	dc[ppgpp_index] = v_rela + v_spot_syn - v_spot_deg - mu * ppgpp
	dc[rnap_bound_index] = v_rnap_binding - v_rnap_completion - mu * rnap_bound
	dc[biomass_index] = mu

	# import ipdb; ipdb.set_trace()
	return dc

# initial conditions - TODO - initialize
co = np.zeros(n_molecules)
co[aa_indices] = 100
co[trna_indices] = 5
co[trna_charged_indices] = 4.5
co[mrna_indices] = 0.1
co[protein_indices] = [0.05, 0.01, 5, 15, 1, 1]  # (RelA, SpoT, RNAP, ribosome, metabolism, synthetases)
co[ppgpp_index] = 10
co[rnap_bound_index] = 1
co[biomass_index] = 1000
t = range(1000)

# solve
sol = odeint(dcdt, co, t)

# plot solution
n_subplots = 4
plt.subplot(n_subplots,1,1)
plt.plot(t, sol[:,ppgpp_index])
plt.ylabel('[ppGpp]')

plt.subplot(n_subplots,1,2)
plt.plot(t, sol[:,protein_indices])
plt.ylabel('[protein]')
plt.legend(['RelA', 'SpoT', 'RNAP', 'ribosome', 'metabolism', 'synthetases'], fontsize=4)

plt.subplot(n_subplots,1,3)
plt.plot(t, sol[:,aa_indices])
plt.ylabel('[AA]')

plt.subplot(n_subplots,1,4)
plt.plot(t, sol[:,trna_charged_indices])
plt.ylabel('[charged tRNA]')
plt.show()

# import ipdb; ipdb.set_trace()