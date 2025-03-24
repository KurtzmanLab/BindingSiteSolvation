# BindingSiteSolvation

This repository provides input files and results for the study of **Solvation Thermodynamic Costs of Cognate Binding Site Formation**. The content is organized into the following folders:


#### **Data**
This folder contains:
- GIST integration results
- SSTMap hydration site analysis results
of 34 DUD-E systems.

#### **MD_Input_files**
This folder contains:
- MD input files for both *rigid* and *flexible* simulations.
- `.prmtop` and `.rst7` files for each system.

  File Details:
  - *Rigid target inputs:*
    - `target_min.prmtop`, `target_min.rst7`
  - *Flexible target inputs:*
    - `target.prmtop`, `target.rst7`

#### **GIST_Input_files**
This folder contains:
- GIST input files for 34 systems.
- The topology and trajectory files for rigid and flexible systems were used as inputs, respectively.

#### **Site_analysis_tools**
This folder contains:
- Python code to analyze binding sites.
- Number of hydrogen bonds, Number of water neighbors, Conformation analysis of residues.
