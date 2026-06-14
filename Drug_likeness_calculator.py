import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, PandasTools
import pandas as pd

# Helper function to check Ghose filter
def check_ghose(mol, mw, logp, mr, num_atoms):
    return (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (40 <= mr <= 130) and (20 <= num_atoms <= 70)

# Helper function to check Veber rule
def check_veber(rot_bonds, tpsa):
    return (rot_bonds <= 10) and (tpsa <= 140)

# Helper function to check Egan rule
def check_egan(tpsa, logp):
    return (tpsa <= 132) and (logp <= 6)

# Helper function to check Muegge filter
def check_muegge(mol, mw, logp, hba, hbd, rot_bonds, tpsa):
    # Allowed atomic numbers: C(6), H(1), O(8), N(7), S(16), P(15), F(9), Cl(17), Br(35), I(53)
    allowed_atoms = {6, 1, 8, 7, 16, 15, 9, 17, 35, 53}
    atom_check = all(atom.GetAtomicNum() in allowed_atoms for atom in mol.GetAtoms())
    ring_count = Lipinski.RingCount(mol)
    return (200 <= mw <= 600) and (-2 <= logp <= 5) and (hba <= 10) and (hbd <= 5) and (rot_bonds <= 15) and (ring_count <= 4) and atom_check and (tpsa <= 150)

# Helper function to check Rule of 3
def check_rule_of_3(mw, logp, hbd, hba, rot_bonds):
    return (mw < 300) and (logp <= 3) and (hbd <= 3) and (hba <= 3) and (rot_bonds <= 3)

def main():
    if len(sys.argv) != 3:
        print("Usage: python Drug_likeness.py input.smi output.csv")
        return
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Read SMILES file
    df = pd.read_csv(input_file, sep='\t', header=None, names=['SMILES', 'Name'])
    
    # Calculate properties
    results = []
    for smi, name in zip(df['SMILES'], df['Name']):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                # Basic Lipinski descriptors
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                
                # Additional descriptors
                heavy_atoms = Lipinski.HeavyAtomCount(mol)
                rot_bonds = Lipinski.NumRotatableBonds(mol)
                tpsa = Descriptors.TPSA(mol)
                mr = Descriptors.MolMR(mol)
                num_atoms = mol.GetNumAtoms()
                ring_count = Lipinski.RingCount(mol)
                
                # Check Lipinski rules
                rule1 = mw <= 500
                rule2 = logp <= 5
                rule3 = hbd <= 5
                rule4 = hba <= 10
                passes_lipinski = rule1 and rule2 and rule3 and rule4
                
                # Check other filters
                passes_ghose = check_ghose(mol, mw, logp, mr, num_atoms)
                passes_veber = check_veber(rot_bonds, tpsa)
                passes_egan = check_egan(tpsa, logp)
                passes_muegge = check_muegge(mol, mw, logp, hba, hbd, rot_bonds, tpsa)
                passes_rule3 = check_rule_of_3(mw, logp, hbd, hba, rot_bonds)
                
                results.append({
                    'Name': name,
                    'SMILES': smi,
                    'MW': round(mw, 2),
                    'LogP': round(logp, 2),
                    'HBD': hbd,
                    'HBA': hba,
                    'Heavy atoms': heavy_atoms,
                    'Num of atoms': num_atoms,
                    'Rotatable bonds': rot_bonds,
                    'Total Polar Surface area': tpsa,
                    'Molar Refraction': mr,
                    'Passes_Lipinski': passes_lipinski,
                    'Ring Count': ring_count,
                    'Passes_Ghose': passes_ghose,
                    'Passes_Veber': passes_veber,
                    'Passes_Egan': passes_egan,
                    'Passes_Muegge': passes_muegge,
                    'Passes_RuleOf3': passes_rule3
                })
        except:
            continue
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save to CSV
    results_df.to_csv(output_file, index=False)
    print(f"Processed {len(results_df)} molecules. Results saved to {output_file}")

if __name__ == "__main__":
    main()