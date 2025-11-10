import os, sys
import json
from Bio import SeqIO
import argparse
import subprocess
import re
from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from collections import defaultdict
sys.path.append('/data/alexandre/scripts/colabfold_updated/running_predictions/alphafold3')
from utils_af3 import (
    # ACCEPT_DEFAULT_TERMS,
    # DEFAULT_API_SERVER,
    # NO_GPU_FOUND,
    # CIF_REVISION_DATE,
    # get_commit,
    # safe_filename,
    # setup_logging,
    CFMMCIFIO
    )

def AAContentPDB(pdb):
    aa31 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    chain_sequences = {}

    with open(pdb, 'r') as filein:
        for line in filein:
            if line.startswith('ATOM') and 'CA' in line and (line[16] == 'A' or line[16] == ' '):
                chain_id = line[21]
                if chain_id not in chain_sequences:
                    chain_sequences[chain_id] = ''
                res = line[17:20]
                aa = aa31[res]
                chain_sequences[chain_id] += aa

    return chain_sequences

def read_fasta(file_path, chain_id):
    with open(file_path, 'r') as file:
        fasta_content = file.read().strip()
    sequences = {}
    # chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    chain_idx=0
    entries = fasta_content.split('>')[1:]
    # getting chain sequence
    first_entry=entries[0]
    elements = first_entry.split('\n')
    sequence = elements[1].strip()
    for chain_sequence in sequence.split(':'):
        # # getting chain ID
        # chain_id = chain_ids[chain_idx]
        # relating chain ID to seq
        sequences[chain_id] = chain_sequence
        chain_idx+=1
    return sequences

# Copied from batch.py
modified_mapping = {
  "MSE" : "MET", "MLY" : "LYS", "FME" : "MET", "HYP" : "PRO",
  "TPO" : "THR", "CSO" : "CYS", "SEP" : "SER", "M3L" : "LYS",
  "HSK" : "HIS", "SAC" : "SER", "PCA" : "GLU", "DAL" : "ALA",
  "CME" : "CYS", "CSD" : "CYS", "OCS" : "CYS", "DPR" : "PRO",
  "B3K" : "LYS", "ALY" : "LYS", "YCM" : "CYS", "MLZ" : "LYS",
  "4BF" : "TYR", "KCX" : "LYS", "B3E" : "GLU", "B3D" : "ASP",
  "HZP" : "PRO", "CSX" : "CYS", "BAL" : "ALA", "HIC" : "HIS",
  "DBZ" : "ALA", "DCY" : "CYS", "DVA" : "VAL", "NLE" : "LEU",
  "SMC" : "CYS", "AGM" : "ARG", "B3A" : "ALA", "DAS" : "ASP",
  "DLY" : "LYS", "DSN" : "SER", "DTH" : "THR", "GL3" : "GLY",
  "HY3" : "PRO", "LLP" : "LYS", "MGN" : "GLN", "MHS" : "HIS",
  "TRQ" : "TRP", "B3Y" : "TYR", "PHI" : "PHE", "PTR" : "TYR",
  "TYS" : "TYR", "IAS" : "ASP", "GPL" : "LYS", "KYN" : "TRP",
  "CSD" : "CYS", "SEC" : "CYS"
}

# Copied from batch.py
class ReplaceOrRemoveHetatmSelect(Select):
  def accept_residue(self, residue):
    hetfield, _, _ = residue.get_id()
    if hetfield != " ":
      if residue.resname in modified_mapping:
        # set unmodified resname
        residue.resname = modified_mapping[residue.resname]
        # clear hetatm flag
        residue._id = (" ", residue._id[1], " ")
        t = residue.full_id
        residue.full_id = (t[0], t[1], t[2], residue._id)
        return 1
      return 0
    else:
      return 1

# Adapted from batch.py
def convert_pdb_to_mmcif(pdb_file):
    """Adapted from batch.py - convert existing pdb files into mmcif with the required poly_seq and revision_date"""
    pdb_dir = os.path.dirname(pdb_file)
    pdb_basename = os.path.basename(pdb_file)
    pdb_name = os.path.splitext(pdb_basename)[0]
    # parse PDB
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(pdb_name, pdb_file)
    # convert PDB to CIF and save it
    cif_file = os.path.join(pdb_dir, f"{pdb_name}.cif")
    cif_io = CFMMCIFIO()
    cif_io.set_structure(structure)
    cif_io.save(str(cif_file), ReplaceOrRemoveHetatmSelect())
    return cif_file

def generate_protein_json(input_name, output_dir, a3m_files, templates_path, model_seeds, msa_type='expert'):
    """
    Generate a JSON file for the protein structure input based on the specified requirements.

    Parameters:
    - input_name (str): input pdbfile
    - output_file (str): Path to save the JSON file.
    - fasta_file (str): Path to the input FASTA file containing chain sequences.
    - a3m_paths (list): List of paths to the unpaired MSA files (A,B,C...).
    - cif_path (str): Path to the CIF template file.
    - model_seeds (list): List of model seeds (e.g., [1, 2, 3, 4, 5]).
    """

    # Build the initial JSON structure
    structure = {
        "name":input_name,
        "sequences": [],
        "modelSeeds": list(range(1, int(model_seeds[0]) + 1)),
        "dialect": "alphafold3",
        "version": 1
    }

    # Converting template from PDB to mmCIF if PDB is provided
    cifs_path_dict = {}
    for template_path in templates_path:
        template_chain = os.path.basename(template_path).split('.')[0].split('_')[-1] # from <path>/<pdb_name>_<chain>.pdb (e.g. templates/5kbn_A.pdb) to <chain> (e.g A)
        if template_path.endswith('.pdb'):
            cif_path = convert_pdb_to_mmcif(template_path)
        elif template_path.endswith('.cif'):
            cif_path = template_path
        else:
            exit('Template provided is in unrecognized format. Please provide a PDB or mmCIF file and try again. Stopping script.')
        cifs_path_dict[template_chain]=cif_path

    # Parse the A3M files to extract the input sequences
    sequences = {}

    # Create dictionary for A3M files (chain ID -> a3m files)
    a3m_files_dict = defaultdict(dict)
    for a3m_file in a3m_files:
        filename = os.path.basename(a3m_file)
        match = re.match(r".*_(paired|unpaired)_([A-Z])\.a3m$", filename)
        if not match:
            raise ValueError(f"Unexpected A3M filename format: {filename}")
        msa_category, chain_id = match.groups()
        a3m_files_dict[chain_id][msa_category] = a3m_file
    print(a3m_files_dict)
    print(msa_type)

    # Iterate through each chain and corresponding a3m files
    for chain_id, current_a3m_files in a3m_files_dict.items():
        num_paired_a3m_files = int('paired' in current_a3m_files)
        num_unpaired_a3m_files = int('unpaired' in current_a3m_files)
        print (chain_id, num_paired_a3m_files, num_unpaired_a3m_files)

        paired_a3m_file = ""
        paired_msa_key = None

        if msa_type == 'expert':
            if num_paired_a3m_files != 1 or num_unpaired_a3m_files != 1:
                sys.exit(f"Inadequate input a3m files for expert configuration. Exactly 1 paired and 1 unpaired a3m files should have been provided for chain {chain_id}")
            paired_a3m_file = current_a3m_files['paired']
            paired_msa_key = "pairedMsaPath"
            print (paired_msa_key, paired_a3m_file)

        elif msa_type == 'af3_recommended':
            if num_paired_a3m_files != 0 or num_unpaired_a3m_files == 0:
                sys.exit(f"Inadequate input a3m files for AF3 recommended configuration — no paired a3m files expected and at least 1 unpaired a3m file should have been provided for chain {chain_id}")
            paired_a3m_file = ""
            paired_msa_key = "pairedMsa"

        unpaired_a3m_file = current_a3m_files['unpaired']

        # query sequence
        current_sequences = read_fasta(unpaired_a3m_file, chain_id)
        current_sequence = current_sequences[chain_id]

        sequences.update(current_sequences)

        # templates
        if chain_id == 'A':
            templates=[]
        else:
            templateIndices = [i for i in range(len(current_sequence))]
            queryIndices = [i for i in range(len(current_sequence))]
            templates = [
                {
                    "mmcifPath": cifs_path_dict[chain_id],
                    "queryIndices": "[" + ", ".join(map(str, queryIndices)) + "]",
                    "templateIndices": "[" + ", ".join(map(str, templateIndices)) + "]"
                }
            ]

        # Defining current chain entry in JSON input file
        chain_entry = {
            "protein": {
                "id": chain_id,
                "sequence": current_sequence,
                "unpairedMsaPath": unpaired_a3m_file,
                paired_msa_key: paired_a3m_file,
                "templates": templates
            }
        }

        # Appending chain entry to JSON structure
        structure["sequences"].append(chain_entry)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pdbname = input_name[:-4]
    output_file = f"{pdbname}_input_af3.json"
    output_filepath = os.path.join(output_dir, output_file)
    # Write to JSON file with standard indent (4 spaces)
    with open(output_filepath, "w") as f:
        json.dump(structure, f, indent=4, separators=(',', ': '))

    # "Hack" - remove outter commas from the string-printed lists
    subprocess.run(['sed', '-i', 's/\"\[/\[/g', output_filepath])
    subprocess.run(['sed', '-i', 's/\]\"/\]/g', output_filepath])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a JSON file for protein structure input.")
    parser.add_argument("-i", "--input_name", type=str, default='test', help="input pdbfile")
    parser.add_argument("-o", "--output_dir", type=str, help="Path to save the JSON file.")
    parser.add_argument("-a3m_files", "--a3m_files", nargs='+', required=True, type=str, help="Path to the unpaired MSA file.")
    parser.add_argument("-msa_type", "--msa_type", type=str, help="Type of MSA provided.")
    parser.add_argument("-templates", "--templates_path", nargs='+', type=str, help="Path to the CIF template file.")
    parser.add_argument("-seeds", "--model_seeds", nargs='+', type=int, default=[1, 2, 3, 4, 5], help="List of model seeds (e.g., 1 2 3 4 5).")

    args = parser.parse_args()

    generate_protein_json(args.input_name, args.output_dir, args.a3m_files, args.templates_path, args.model_seeds, args.msa_type)
