File sources:
- gene_ids.tsv: validation/ecoli/flat/geneIDs.tsv
- ppgpp_fc.tsv: supplemental data (Dataset_S01) from Sanchez-Vazquez et al. Genome-wide effects on Escherichia coli transcription from ppGpp binding to its two sites on RNA polymerase. PNAS. 116(17) 8310-8319. 2019. https://www.pnas.org/highwire/filestream/859108/field_highwire_adjunct_files/1/pnas.1819682116.sd01.xlsx
- ppgpp_regulation.tsv: manually curated from https://ecocyc.org/compound?orgid=ECOLI&id=GUANOSINE-5DP-3DP#tab=REG, to be added to reconstruction/ecoli/flat/
- sim_data.cp: simData.cPickle output from `python runscripts/manual/runParca.py -c8` on `191002-ppgpp-reg 9b1bfea4f`