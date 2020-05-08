Explore data source of foldChanges.tsv in reconstruction.

File sources:
- kb.tsv: output/data EcoliFC/KB.txt (from EcoliFoldChanges repo - 8580aae)
- shifts.tsv: data/EcoMAC/experimentalShift.tsv (from EcoliFoldChanges repo - 8580aae)
- gene_names.tsv: data/EcoMAC/geneName.tsv (from EcoliFoldChanges repo - 8580aae)
- fold_changes.tsv: reconstruction/ecoli/flat/foldChanges.tsv (from wcEcoli repo - 143cf74856)
- fc_single_shift.tsv: output/fc_single_shift.tsv (from running src/condition_graph.py in FoldChangeExplorer repo - 7824581)
- fc_all_shifts.tsv: output/fc_all_shifts.tsv (from running src/condition_graph.py in FoldChangeExplorer repo - 7824581)
- fc_aggregated_shifts.tsv: output/fc_aggregated_shifts.tsv (from running src/condition_graph.py in FoldChangeExplorer repo - 7824581)

Running:
- `./calculate_fc.py`: attempts to use best information available (regulation direction discrepancies)
- `./calculate_fc.py -m`: attempts to match expected processing of source data (relies on unknown method of determining regulation direction fater processing and compares based on assuming only positive direction)

Notes:
Original data appears to ignore direction of regulation in data analysis.
Mean and std are calculated based on abs of fold change and direction is
probably set afterwards based on curated regulation.  This can lead to issues
with ambiguous regulation where annotation point to positive and negative
regulation (eg ArcA -> sdhCDAB is positive even though fold change data points
to a negative regulation).

>90% of data appears to be consistent (82/938 fold changes do not match with -m option and ignoring fold changes missing in wcm).
