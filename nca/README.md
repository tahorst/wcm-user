Exploration of applying network component analysis (NCA) to expression data to get fold change data for the model.

Relevant papers:
- Liao et al. Network component analysis: Reconstruction of regulatory signals in biological systems. PNAS. 2003.
- Chang et al. Fast network component analysis (FastNCA) for gene regulatory network reconstruction from microarray data. Bioinformatics. 2008.
- Noor et al. ROBNCA: robust network component analysis for recovering transcription factor activities. Bioinformatics. 2013.
- Jayavelu et al. Iterative sub-network component analysis enables reconstruction of large scale genetic networks. BMC Bioinformatics. 2015.
- Carrera et al. An integrative, multi‐scale, genome‐wide model reveals the phenotypic landscape of Escherichia coli. MSB. 2014.

Setup (optional depending on desired regulatory network topology):
Run `./download_regulondb.sh` to download data from regulonDB that is not committed.

TODO:
- Normalize sequencing datasets (might need Lowess normalization for each sample and quantile normalization between all samples) to use datasets in addition to EcoMAC.tsv
- Consider new fold change calculation for linearized input data (current calculation assumes input is on a log scale and might not be appropriate for linear predictions)
