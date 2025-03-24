## Setup and Run

### 1. Create Environment

First, create a new conda environment using the `environment.yml` file provided in the project. Navigate to the project folder and run the following command:

```
conda env create -f environment.yml
```
This will set up the necessary environment for the project.


### 2. Your Specific File Paths and Names
In the `scripts/run.py` file, you'll need to specify the following:

- **Target file:** Specify the path to the target protein file (in PDB format).
- **Residue indices:** Update the indices for the residues (as per your PyMOL representation).
- **Trajectory frames to analyze:** Specify the frames you want to analyze from your trajectory file.
- **Files you need:** Rigid & flexible time-average files (.pdb), Topology (.ptmtop) & trajectory (.nc) files

To modify the script, search for the following sections in run.py and make adjustments to suit your files.


### 3. Navigate to the Project Directory
Next, navigate to the project directory:
```
cd path_to_project/Site_analysis_tools
```

### 4. Run the Script
After making the necessary modifications to run.py, you can execute the script by running:
```
python -m scripts.run
```

### 5. Results
After the script finishes running, the results will be saved in the `results` folder within the project directory.
