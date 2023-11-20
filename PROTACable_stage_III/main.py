import os, sys
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from elabmol import ElabMols
import pymol
from multiprocessing import Pool, cpu_count, Manager
from openbabel import pybel
ob = pybel.ob
from openbabel import openbabel

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process pdb files.')
    parser.add_argument('POCKET', type=str, help='Path to the POCKET pdb file.')
    parser.add_argument('POI_LIGAND', type=str, help='Path to the POI LIGAND pdb file.')
    parser.add_argument('E3_LIGAND', type=str, help='Path to the E3 LIGAND pdb file.')
    parser.add_argument('--LINKER', type=str, default=None, help='Path to the LINKER pdb file.')
    return parser.parse_args()

def monitor_termination_flag(shared_dict, pool):
    while True:
        if shared_dict.get('terminate'):
            pool.terminate()
            break

def correct_valency_issues(molecule):
    typical_valencies = {
        'C': [4], 'N': [3], 'O': [2], 'H': [1], 
        'P': [3, 5], 'S': [2, 4, 6], 'B': [3]
    }

    corrected_molecule = Chem.RWMol(molecule)

    atoms_to_remove = []

    for atom in corrected_molecule.GetAtoms():
        element = atom.GetSymbol()
        if element in typical_valencies:
            if atom.GetDegree() not in typical_valencies[element]:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H':
                        atoms_to_remove.append(neighbor.GetIdx())

    for idx in sorted(atoms_to_remove, reverse=True):
        corrected_molecule.RemoveAtom(idx)

    return corrected_molecule.GetMol()

def merge_with_pymol(pocket_pdb_path, molecule_pdb_path, output_path):
    pymol.finish_launching(['pymol', '-qc'])
    pymol.cmd.delete("all")

    pymol.cmd.load(pocket_pdb_path, 'pocket')
    pymol.cmd.load(molecule_pdb_path, 'molecule')

    pymol.cmd.save(output_path, 'all')

    # pymol.cmd.quit()

def convert_pdb_to_mol2(pdb_input_path, mol2_output_path=None):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, pdb_input_path)
    mol.AddHydrogens()

    if not mol2_output_path:
        mol2_output_path = f"{os.path.splitext(pdb_input_path)[0]}.mol2"

    obConversion.WriteFile(mol, mol2_output_path)
    return mol2_output_path

def process_linker(linker_full_path, pocket, poi_ligand, e3_ligand, output_dir, shared_dict):
    em = ElabMols()
    result_mol, skip = em.elaborate_two_mols(pocket, poi_ligand, e3_ligand, frag_path=linker_full_path)


    if skip:
        print(f"Pocket confirmation {os.path.basename(pocket)} was rejected")
        terminate_event.set()
        return

    if result_mol is None:
        print(f"Linker {os.path.basename(linker_full_path)} was rejected")
        return

   # print("writing")
    result_mol = correct_valency_issues(result_mol)
    temp_output_path = os.path.join(output_dir, f"temp_{os.path.splitext(os.path.basename(poi_ligand))[0]}_{os.path.splitext(os.path.basename(e3_ligand))[0]}_{os.path.splitext(os.path.basename(linker_full_path))[0]}.pdb")
    Chem.MolToPDBFile(result_mol, temp_output_path)

    output_name_base = f"{os.path.splitext(os.path.basename(poi_ligand))[0]}_{os.path.splitext(os.path.basename(e3_ligand))[0]}_{os.path.splitext(os.path.basename(linker_full_path))[0]}"
    output_path_pdb = os.path.join(output_dir, output_name_base + ".pdb")

    merge_with_pymol(pocket, temp_output_path, output_path_pdb)

    ligands_dir = os.path.join(output_dir, 'ligands')
    if not os.path.exists(ligands_dir):
        os.makedirs(ligands_dir)

    output_path_mol2 = os.path.join(ligands_dir, output_name_base + ".mol2")
    output_path_sdf = os.path.join(ligands_dir, output_name_base + ".sdf")

    convert_pdb_to_mol2(temp_output_path, output_path_mol2)
    AllChem.MolToMolFile(result_mol, output_path_sdf)

    os.remove(temp_output_path)
    print(f"Linker {os.path.basename(linker_full_path)} was accepted. Results saved in {output_dir}")

def main():
    args = parse_arguments()

    output_dir = "output_stage_3_4"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"Output directory created at: {os.path.abspath(output_dir)}")

    manager = Manager()
    terminate_event = manager.Event()

    linkers_to_process = []
    if args.LINKER:
        linkers_to_process.append(args.LINKER)
    else:
        for linker_path in os.listdir(f"{os.environ['PROTACable']}/PROTACable_stage_3/linker_library/filtered_pdbs/"):
            if linker_path.endswith('.pdb'):
                linkers_to_process.append(os.path.join(f"{os.environ['PROTACable']}/PROTACable_stage_3/linker_library/filtered_pdbs/", linker_path))

    process_args = [(linker, args.POCKET, args.POI_LIGAND, args.E3_LIGAND, output_dir, terminate_event) for linker in linkers_to_process]

    with Pool(processes=24) as pool:
        async_result = pool.starmap_async(process_linker, tqdm(process_args, desc="Processing linkers"))
        while not async_result.ready():
            if terminate_event.is_set():
                pool.terminate()
                break
            async_result.wait(1)  # Wait for 1 second before checking again

    if terminate_event.is_set():
        print("Termination flag set. Exiting...")
        sys.exit(1)

if __name__ == "__main__":
    main()
