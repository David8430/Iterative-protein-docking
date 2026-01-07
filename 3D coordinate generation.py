from os.path import isfile, splitext
from fontTools.misc.cython import returns
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from os import path, listdir
from rdkit.Chem.AllChem import EmbedMolecule
from rdkit.Chem.AllChem import AssignStereochemistry
from rdkit.Chem.AllChem import UFFOptimizeMolecule
from rdkit.Chem.EnumerateStereoisomers import GetStereoisomerCount
from rdkit.VLib.NodeLib.demo import output


def generate_3d_sdf(molecule_b):
    molecule_b = Chem.AddHs(molecule_b, addCoords = True)

    params = AllChem.EmbedParameters()
    params.randomSeed = 5
    params.useRandomCoords = True
    Chem.AllChem.EmbedMolecule(molecule_b, params)

    if molecule_b.GetNumConformers() == 0:
        raise ValueError('3D embedding failed')

    AllChem.UFFOptimizeMolecule(molecule_b)

    output_name = input('What name do you want to save it under? ') + '.sdf'
    script_location = os.path.dirname(os.path.abspath(__file__))
    if 'StructureOutput' not in listdir(script_location):
        os.mkdir(os.path.join(script_location, 'StructureOutput'))
    writer_tool = Chem.SDWriter(os.path.join(script_location, 'StructureOutput', output_name))
    writer_tool.write(molecule_b)
    writer_tool.close()

def smiles_input():
    mol_smiles = input("Please enter the SMILES string: ")
    molecule_structure = Chem.MolFromSmiles(mol_smiles, sanitize=True)

    try:
        Chem.SanitizeMol(molecule_structure)
    except ValueError:
        print('Invalid molecule input.')
    else:
        if molecule_structure is None:
            print("invalid input")
            return None
        else:
            print("SMILES accepted")
            if GetStereoisomerCount(molecule_structure) == 1:
                return molecule_structure
            else:
                print(f'Isomer not precisely defined.'
                      f' A total of {GetStereoisomerCount(molecule_structure)} possibilities were found.'
                      f' Please remove ambiguity from input.')
                return None

def file_input():
    mol_path = input("Please enter the file's full path: ")
    if os.path.isfile(mol_path) and splitext(mol_path)[1] == '.sdf':
        molecule_structure = Chem.MolFromMolFile(mol_path, sanitize=True, removeHs=False)

        try:
            Chem.SanitizeMol(molecule_structure)
        except ValueError:
            print('Invalid molecule input.')
        else:
            if molecule_structure is not None:
                if GetStereoisomerCount(molecule_structure) == 1:
                    return molecule_structure
                else:
                    print(f'Isomer not precisely defined.'
                          f' A total of {GetStereoisomerCount(molecule_structure)} possibilities were found.'
                          f' Please remove ambiguity from input.')
                    return None
            else:
                print('Molecule creation failed.')
                return None
    else:
        print('Incorrect/non-existent input')
        return None

def inchi_input():
    mol_inchi = input("Please enter the InChI string: ")
    molecule_structure = Chem.inchi.MolFromInchi(mol_inchi, sanitize=True, removeHs=False)

    try:
        Chem.SanitizeMol(molecule_structure)
    except ValueError:
        print('Invalid molecule input.')
    else:
        if molecule_structure is None:
            print("invalid input")
            return None
        else:
            print("InChI accepted")
            if GetStereoisomerCount(molecule_structure) == 1:
                return molecule_structure
            else:
                print(f'Isomer not precisely defined.'
                      f' A total of {GetStereoisomerCount(molecule_structure)} possibilities were found.'
                      f' Please remove ambiguity from input.')
                return None

def take_input():
    while True:
        input_mode = input("Please choose from the following input possibilities:\n1-SMILES string\n2-.sdf file\n3-InChI code\n")
        if input_mode == '1':
            target_molecule = smiles_input()
            if target_molecule is not None:
                return target_molecule
        elif input_mode == '2':
            target_molecule = file_input()
            if target_molecule is not None:
                return target_molecule
        elif input_mode == '3':
            target_molecule = inchi_input()
            if target_molecule is not None:
                return target_molecule
        else:
            print("Invalid input choice.")

# Run Script
if __name__ == '__main__':
    print('Welcome to the 3D molecular structure generator.')
    while True:
        processed_molecule = take_input()

        generate_3d_sdf(processed_molecule)

        if input("Do you wish to continue? (y/n)") == 'n':
            break