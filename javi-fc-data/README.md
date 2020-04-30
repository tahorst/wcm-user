Explore data source of foldChanges.tsv in reconstruction.

File sources:
- kb.tsv: output/data EcoliFC/KB.txt (from EcoliFoldChanges repo - 8580aae)
- fold_changes.tsv: reconstruction/ecoli/flat/foldChanges.tsv (from wcEcoli repo - 143cf74856)

Notes:
Original data appears to ignore direction of regulation in data analysis.
Mean and std are calculated based on abs of fold change and direction is
probably set afterwards based on curated regulation.  This can lead to issues
with ambiguous regulation where annotation point to positive and negative
regulation (eg ArcA -> sdhCDAB is positive even though fold change data points
to a negative regulation).

>90% of data appears to be consistent (82/938 fold changes do not match, ignoring fold changes missing in wcm).
