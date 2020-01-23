Scripts to explore adding additional kinetics data to metabolism.  Currently have 1268 constraints annotated but only use 431 of them (some of which might be split in two from a single spreadsheet constraint if multiple kcats exist for different metabolites).

Requires:
- data/raw_data.cp
- data/sim_data.cp

Tested with raw_data.cp and sim_data.cp from parca output files in out/manual/kb/ (rawData.cPickle and simData.cPickle) on range-kinetics branch (f5f512a27) branched from edd286050e with the addition of metabolism_kinetics.tsv flat file.
