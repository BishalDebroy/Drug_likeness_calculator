# Drug Likeness Calculator

A Python tool that evaluates multiple drug‑likeness rules for a set of molecules given as SMILES strings.  
It uses **RDKit** to compute molecular descriptors and checks compliance with the following filters:

| Filter              | Criteria                                                                                                 |
|---------------------|----------------------------------------------------------------------------------------------------------|
| Lipinski (Ro5)      | MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10                                                                   |
| Ghose               | MW 160–480, LogP –0.4 to 5.6, Molar refractivity 40–130, total atoms 20–70                              |
| Veber               | Rotatable bonds ≤ 10, TPSA ≤ 140 Å²                                                                     |
| Egan                | TPSA ≤ 132 Å², LogP ≤ 6                                                                                 |
| Muegge              | MW 200–600, LogP –2 to 5, HBA ≤10, HBD ≤5, rotatable bonds ≤15, rings ≤4, allowed atoms, TPSA ≤150 Å² |
| Rule of 3 (fragments) | MW < 300, LogP ≤3, HBD ≤3, HBA ≤3, rotatable bonds ≤3                                                 |

# Requirements
1. Linux or WSL
2. rdkit>=2023.03.1
3. pandas>=1.5.0

# Suggestion
Create an rdkit environment with python. Then install the requirements
```
conda create -n rdkit_env python==3.13
conda activate rdkit_env
pip install rdkit pandas
conda deactivate
```


## Input format
Tab‑separated file (`.smi` or `.txt`) with **two columns** and **no header**: Check examples


# How to run the program
1. Clone the repository ina suitable working directory in your linux
   ```
   git clone https://github.com/BishalDebroy/Drug_likeness_calculator.git
   ```
   
2. Go into the directory and replace the content of `SMILE.smi` with your list of molecule(s)
3. execute the program
```
conda activate rdkit_env
python drug_likeness_calculator.py SMILES.smi dlc_output.csv
```

# To update the installed packages
```
conda clean --all
conda upgrade -n rdkit_env --all
```

# To update conda
```
conda update -n base conda
```
