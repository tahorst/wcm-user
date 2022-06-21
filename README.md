# wcm-user
A collection of ad hoc and exploratory scripts and side projects related to the [_E. coli_ whole-cell model](https://github.com/CovertLab/WholeCellEcoliRelease).
This repo was included as the `user/` sub-directory in the whole-cell model repo during development throughout my PhD and depends on the environment setup and relative paths in that repo.
Some large data files are not committed so not all scripts are guaranteed to run without regenerating the data files.
Additionally, some scripts are dependent on older versions of the computing environment (Python2 and previous dependencies)
and have not been updated as the environment has been updated, which may lead to issues if attempting to run.
Some subdirectories include a README to give an overview of the files included.

This collection of code spans many years of development while I learned better coding practices along the way.
Many of the newer subdirectories will be more representative of my current approach to coding such as:
- `grc/` - exploration into growth rate control
- `growth-paper/` - analysis scripts for a manuscript under review
- `metabolism-debugger/` - framework to make debugging the metabolism submodel in simulations easier
- `nca/` - applying network component analysis (NCA) to gene expression data
