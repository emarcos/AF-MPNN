#!/usr/bin/env python3
from prody import *
from argparse import *
import os, sys
# from pyrosetta import init, pose_from_pdb
# #from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
# from pyrosetta import *
# from pyrosetta.rosetta import *
# from pyrosetta.rosetta.core.select.residue_selector import *
from argparse import ArgumentParser
import re

def process_mpnn_output(input_file, output_file, chain, fasta, verbose=False):

    """
    Process the ProteinMPNN output file to
    update sequence titles.
    """

    if verbose:
        print(f"Processing ProteinMPNN paired MSA file: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")

    nseq = 0
    with open(input_file, 'r') as f, open(output_file, 'w') as o:
        for line_number, line in enumerate(f):
            line = line.strip()
            if line_number < 1:
                #o.write(line[:(len(fasta)])
                continue
            if not line.startswith('>'):
                if chain == 'A':
                    line = line.split('/')[0]+'\n'
                elif chain == 'B':
                    line = line.split('/')[1]+'\n'
                elif chain == 'C':
                    line = line.split('/')[2]+'\n'
            else:
                line = f'>seq_{nseq:03d}\n'
                nseq += 1
            o.write(line)


    if verbose:
        print(f"ProteinMPNN paired MSA processing completed.")

def cut_aligned_seq(chains, aligned_seq):

    insertion_pattern = r'^[a-z]$' # lowercase letters

    # dict containing output trimmed seqs
    trim_dict = {chain_id: '' for chain_id in chains}

    # initial params
    count = -1
    initial_limit = 0

    # iterating through initial chains dict
    for chain_id, chain_seq in chains.items():
        
        if not chain_seq:
            trim_dict[chain_id]=None
            continue
            
        print(f'Current chain: {chain_id}')

        limit = len(chain_seq)+initial_limit # calc current length limit
        aligned_seq_cut = ''

        # iterate through aligned seq while we are under current chain's length lim
        while count < limit-1:
            count+=1
            r = aligned_seq[count] # obtaining current residue under current count
            aligned_seq_cut+=r

            print(count)
            print(r)

            # if char is a lowercase letter, we expand the limit
            if re.match(insertion_pattern, r):
                limit+=1
                print(limit)

        # when the limit is reached for the current chain, we add the cut aligned seq to the output dict
        trim_dict[chain_id]=aligned_seq_cut
        # and we set a new initial limit for the next chain if there is one
        initial_limit=limit

    return trim_dict

def write_paired_a3m_files(pred_type, multimer_msa_a3m, 
                          fastaA, a3m_paired_a,
                          fastaB, a3m_paired_b, 
                          fastaC=None, a3m_paired_c=None,
                          verbose=False):
    
    chains_seqs = {'A': fastaA,
            'B': fastaB,
            'C': fastaC}
    
    if args.verbose:
        print(f"Generating paired A3M files for {pred_type} mode.")

    # adding seq titles
    # chain A
    paired_out_a = open(a3m_paired_a, 'w')
    paired_out_a.write(f'>101\n{fastaA}\n')   

    # chain B
    paired_out_b = open(a3m_paired_b, 'w')
    paired_out_b.write(f'>102\n{fastaB}\n')

    if a3m_paired_c:
        paired_out_c = open(a3m_paired_c, 'w')
        paired_out_c.write(f'>103\n{fastaC}\n')

    # adding paired MSA from ProteinMPNN
    if pred_type == 'mpnn_custom_msa':

        mpnnfile = os.path.join(os.getcwd(), f'seqs/{pdbfilename}.fa')
        processed_mpnnfile_A = os.path.join(os.getcwd(), f'seqs/processed_{pdbfilename}_A.fa')
        processed_mpnnfile_B = os.path.join(os.getcwd(), f'seqs/processed_{pdbfilename}_B.fa')
        process_mpnn_output(mpnnfile, processed_mpnnfile_A, 'A', fastaA, verbose=args.verbose)
        process_mpnn_output(mpnnfile, processed_mpnnfile_B, 'B', fastaB, verbose=args.verbose)

        with open(processed_mpnnfile_A, 'r') as mf_a:
            for line_number, line in enumerate(mf_a):
                if line_number >= 1:
                    paired_out_a.write(f'{line}')
        with open(processed_mpnnfile_B, 'r') as mf_b:
            for line_number, line in enumerate(mf_b):
                if line_number >= 1:
                    paired_out_b.write(f'{line}')

        if fastaC: 
            processed_mpnnfile_C = os.path.join(os.getcwd(), f'seqs/processed_{pdbfilename}_C.fa')
            process_mpnn_output(mpnnfile, processed_mpnnfile_C, 'C', fastaC, verbose=args.verbose)

            with open(processed_mpnnfile_C, 'r') as mf_c:
                for line_number, line in enumerate(mf_c):
                    if line_number >= 1:
                        paired_out_c.write(f'{line}')
    
    elif pred_type == 'default':
        with open(multimer_msa_a3m, 'r') as f:
            inside_block_paired=False
            nseq=0
            for l in f:
                l = l.strip()
                if ">101" in l and "102" in l:
                    inside_block_paired = True
                    continue
                if l == ">101":
                    break
                if inside_block_paired:
                    if l.startswith('>'):
                        pro_l = f'>seq_{nseq:06d}\n' # we add e.g. 6 to ensure that we are not running out of padding
                        paired_out_a.write(pro_l)
                        paired_out_b.write(pro_l)
                        paired_out_c.write(pro_l) if a3m_paired_c else None
                        nseq += 1
                    else:
                        if nseq==0:
                            continue
                        trimmed_chains = cut_aligned_seq(chains_seqs, l) 
                        paired_out_a.write(f"{trimmed_chains['A']}\n")
                        paired_out_b.write(f"{trimmed_chains['B']}\n")
                        paired_out_c.write(f"{trimmed_chains['C']}\n") if a3m_paired_c else None

    paired_out_a.close()
    paired_out_b.close()
    paired_out_c.close() if a3m_paired_c else None
                        

def write_a3m_file(chain, fastaA, fastaB, fastaC, multimer_msa_a3m, output_a3m_unpaired, verbose=False):

    """
    Write the unpaired a3m file for the current chain.
    """

    if verbose:
        print(f"Generating a3m file: {output_a3m_unpaired}")

        ## UNPAIRED MSA ##
        # 3. Appending unpaired MSA a3m file content (monobody and receptor)
        inside_block_unpaired=False
        if args.multimer_msa_a3m:
            with open(multimer_msa_a3m, 'r') as f, open(output_a3m_unpaired, 'w') as unpaired_out:
                if chain == 'A' and fastaC:
                    for l in f:
                        # If we find '>101', start writing lines
                        if l.strip() == ">101":
                            #print ("hello")
                            inside_block_unpaired = True
                            unpaired_out.write(l)  # Write the '>101' line
                            continue
                        # If we find '>102', stop writing lines
                        if l.strip() == ">102":# and inside_block_unpaired:
                            inside_block_unpaired = False
                            break  # Exit the loop as we are done copying
                        # Write lines only if inside the '>101' block
                        if inside_block_unpaired:
                        # Use len_B and len_C to trim the line
                            if l.startswith('>'):
                                unpaired_out.write(l)
                            else:
                                trimmed_line = l[:-(len(fastaB) + len(fastaC)+1)].strip()# Trim the line to the desired length
                                unpaired_out.write(trimmed_line + '\n')  # Write trimmed line              

                if chain == 'A' and not fastaC:
                    for l in f:
                        # If we find '>101', start writing lines
                        if l.strip() == ">101":
                            inside_block_unpaired = True
                            unpaired_out.write(l)  # Write the '>101' line
                            continue
                        # If we find '>102', stop writing lines
                        if l.strip() == ">102":# and inside_block_unpaired:
                            inside_block_unpaired = False
                            break  # Exit the loop as we are done copying
                        # Write lines only if inside the '>101' block
                        if inside_block_unpaired:
                        # Use len_B and len_C to trim the line
                            if l.startswith('>'):
                                unpaired_out.write(l)
                            else:
                                trimmed_line = l[:-(len(fastaB)+1)].strip() # Trim the line to the desired length
                                unpaired_out.write(trimmed_line + '\n')  # Write trimmed line

                if chain == 'B' and fastaC:
                    for l in f:
                        # If we find '>102', start writing lines
                        if l.strip() == ">102":
                            inside_block_unpaired = True
                            unpaired_out.write(l)  # Write the '>102' line
                            continue
                        # If we find '>103', stop writing lines
                        if l.strip() == ">103" and inside_block_unpaired:
                            inside_block_unpaired = False
                            break  # Exit the loop as we are done copying
                        # Write lines only if inside the '>102' block
                        if inside_block_unpaired:
                        # Use len_A and len_C to trim the line
                            if l.startswith('>'):
                                unpaired_out.write(l)
                            else:
                                trimmed_line = l[len(fastaA):-(len(fastaC)+1)].strip()  # Trim the line to the desired length
                                unpaired_out.write(trimmed_line + '\n')  # Write trimmed line

                if chain == 'B' and not fastaC:
                    for l in f:
                        # If we find '>102', start writing lines
                        if l.strip() == ">102":
                            inside_block_unpaired = True
                            unpaired_out.write(l)  # Write the '>102' line
                            continue
                        # If we find '>103', stop writing lines
                        if line.strip() == ">103" and inside_block_unpaired:
                            inside_block_unpaired = False
                            break  # Exit the loop as we are done copying
                        # Write lines only if inside the '>102' block
                        if inside_block_unpaired:
                        # Use len_A and len_C to trim the line
                            if l.startswith('>'):
                                unpaired_out.write(l)
                            else:
                                trimmed_line = l[len(fastaA):].strip() # Trim the line to the desired length
                                unpaired_out.write(trimmed_line + '\n')  # Write trimmed line
                if chain == 'C' and fastaC:
                    for l in f:
                        # If we find '>103', start writing lines
                        if l.strip() == ">103":
                            # print ("hello")
                            inside_block_unpaired = True
                            unpaired_out.write(l)  # Write the '>103' line
                            continue
                        # If we find get to the end of the file
                        if l.strip() == None and inside_block_unpaired:
                            inside_block_unpaired = False
                            break  # Exit the loop as we are done copying
                        # Write lines only if inside the '>103' block
                        if inside_block_unpaired:
                        # Use len_A and len_B to trim the line
                            if l.startswith('>'):
                                unpaired_out.write(l)
                            else:
                                trimmed_line = l[(len(fastaA) + len(fastaB)):].strip() # Trim the line to the desired length
                                unpaired_out.write(trimmed_line + '\n')  # Write trimmed line

    # if verbose:
        # print(f"a3m file generation completed, saved in: {output_a3m_file}")

if __name__ == "__main__":

    parser = ArgumentParser(description=f"""Generate a3m file for AFM-ProteinMPNN prediction.
Example of usage: ./a3m_generator.py -pdb <input_pdb> -fasta <input_fasta> -afm_a3m <input_afm_a3m>""")
    parser.add_argument("-pdb", "--pdbfile", required=True, help="Path to the input PDB file.")
    parser.add_argument("-fasta", "--fastafile", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-afm_a3m", "--multimer_msa_a3m", required=False, help="Path to the input default AlphaFold-Multimer MSA a3m file.")
    parser.add_argument("-pred_type", "--pred_type", type=str, required=True, help="Prediction type (AFM or AFM-ProteinMPNN)")
    parser.add_argument("-o", "--output", help="Path to the desired output file. Defaults to current directory.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode.")

    args = parser.parse_args()

    # 1. Analyzing PDB
    pdbfilename = '.'.join(os.path.basename(args.pdbfile).split('.')[:-1])
    structure = parsePDB(args.pdbfile)
    input_chains = sorted(list(set(chain.getChid() for chain in structure.iterChains())))

    if args.verbose:
        print(f"\nRunning a3m_generator.py on {pdbfilename}...\n")

    if len(input_chains) not in  [2,3]:
        sys.exit(f"Input PDB {args.pdbfile} contains {len(input_chains)} and therefore cannot be processed. Please review and try again.")

    # 2. Analyzing FASTA
    with open(args.fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                fasta_all = line

    fastaA = fasta_all.split(':')[0]
    fastaB = fasta_all.split(':')[1].strip()
    fastaC = None
    num_fasta_seqs = len(fasta_all.split(':'))
    if num_fasta_seqs == 3:
        fastaC = fasta_all.split(':')[2].strip()

    # 3. Quality check before proceeding
    if num_fasta_seqs not in [2,3] or len(input_chains) not in [2,3] or num_fasta_seqs != len(input_chains):
        sys.exit(f"""Input FASTA and PDB files must have the same number of chains (2 or 3).
        - Your FASTA ({args.fastafile}) contains {num_fasta_seqs} chains.
        - Your PDB ({args.pdbfile}) contains {len(input_chains)} chains.
        Please review and try again.""")

    # 5. Generating a3m file
    # PAIRED
    output_a3m_paired_file_A = os.path.join(args.output, f"{pdbfilename}_{args.pred_type}_paired_A.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_{args.pred_type}_paired_A.a3m")
    output_a3m_paired_file_B = os.path.join(args.output, f"{pdbfilename}_{args.pred_type}_paired_B.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_{args.pred_type}_paired_B.a3m")
    if fastaC:
        output_a3m_paired_file_C = os.path.join(args.output, f"{pdbfilename}_{args.pred_type}_paired_C.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_{args.pred_type}_paired_C.a3m")
    else:
        output_a3m_paired_file_C = None

    write_paired_a3m_files(args.pred_type, args.multimer_msa_a3m, 
                          fastaA, output_a3m_paired_file_A,
                          fastaB, output_a3m_paired_file_B, 
                          fastaC, output_a3m_paired_file_C, 
                          args.verbose)
    
    # UNPAIRED
    output_a3m_unpaired_file_A = os.path.join(args.output, f"{pdbfilename}_unpaired_A.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_unpaired_A.a3m")
    write_a3m_file('A', fastaA, fastaB, fastaC, args.multimer_msa_a3m, output_a3m_unpaired_file_A, verbose=args.verbose)

    output_a3m_unpaired_file_B = os.path.join(args.output, f"{pdbfilename}_unpaired_B.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_unpaired_B.a3m")
    write_a3m_file('B', fastaA, fastaB, fastaC, args.multimer_msa_a3m, output_a3m_unpaired_file_B, verbose=args.verbose)
    
    if fastaC:
        output_a3m_unpaired_file_C = os.path.join(args.output, f"{pdbfilename}_unpaired_C.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_unpaired_C.a3m")
        write_a3m_file('C', fastaA, fastaB, fastaC, args.multimer_msa_a3m, output_a3m_unpaired_file_C, verbose=args.verbose)

# end of script
