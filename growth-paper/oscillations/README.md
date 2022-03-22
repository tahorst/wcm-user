*.npy files come from remove_aa_inhibition variants (variant 0 - control, variant 4 -leuA) and show oscillations in the Leu concentration in the LeuA mutant

```
DESC="Remove amino acid inhibition large" PARALLEL_PARCA=1 RUN_AGGREGATE_ANALYSIS=0 \
  N_GENS=16 N_INIT_SIMS=16 \
  VARIANT=remove_aa_inhibition FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=7 \
  TRNA_ATTENUATION=1 PPGPP_REGULATION=1 MECHANISTIC_TRANSLATION_SUPPLY=1 MECHANISTIC_AA_TRANSPORT=1 AA_SUPPLY_IN_CHARGING=1 \
  D_PERIOD_DIVISION=1 MECHANISTIC_REPLISOME=0 TIMESTEP_MAX=1 \
  python runscripts/fireworks/fw_queue.py
```
