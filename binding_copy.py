import Bio.PDB #imports
from Bio.PDB import *
import argparse
import os
import glob
import numpy as np
import string
import pandas as pd
import tkinter as tk
from tkinter import filedialog
from Bio.PDB.Polypeptide import three_to_one, protein_letters_3to1
from Bio import BiopythonWarning
import warnings
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio import Align
import csv
from collections import defaultdict



def is_nucleic_acid(residue):
    """checks if a residue is a nucleic acid

    Args:
        residue (residue class object): accessed through pipeline of structure --> model --> chain --> residue, consists of atoms

    Returns:
        boolean: whether the residue is a nucleic acid

    note: A, T, G, etc represent ribonucleic acids, while DA, DT, DG, etc represent deoxyribonucleic acids
    """
    nucleic_acids = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"} #dictionary of nucleic acids
    return residue.get_resname() in nucleic_acids


def is_amino_acid(residue):
    """checks if residue is amino acid

    Args:
        residue (residue class): see is_nucleic_acid

    Returns:
        boolean: whether the residue is an amino acid
    """
    return (residue.get_resname() in standard_aa_names) #standard_aa_names is a dictionary with all amino acid codons (3 letters)


def calc_amino_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues

    Args:
        residue_one (residue class): first residue to compare
        residue_two (residue class): second residue to compare

    Returns:
        distance: direct distance between the CA of the two residues

    note: these two residues must be amino acids, they are the only residue type with the C-alpha attribute
    """

    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord

    return np.sqrt(np.sum(diff_vector * diff_vector))



def calc_amino_nuc_dist(residue_one, residue_two): #residue one is amino, residue 2 is nucleic
    """ calculates corresponding distances between the CA, phosphate, and base of an amino acid and nucleic acid
    Args:
        residue_one (residue class): first residue to compare, must be nucleic or amino acid
        residue_two (residue class): second residue to compare, must be nucleic or amino acid

    one of the residues must be a nucleic acid and the other an amino acid

    Returns:
        one_p : distance between CA of amino acid and phosphate of nucleic acid
        one_b: distance betweeen CA of amino acid and base of nucleic acid
    """

    nucleic_base_atoms = {
        "A": "N9", "G": "N9", "C": "N1", "T": "N1", "U": "N1",
        "DA": "N9", "DG": "N9", "DC": "N1", "DT": "N1", "DU": "N1" #dictionary with the  appropriate correpsonding molecule based on if the nucleic acid is pyramidine or purine
    }

    DNA_bases = {"DA", "DT", "DG", "DC", "DU"}


    if is_amino_acid(residue_one): #if residue one is the amino acid
        res_one_coord = residue_one["CA"].coord

        if residue_two.get_resname() in DNA_bases:
            if "P" in residue_two:
                nuc_p = residue_two['P'].coord
            else:
                nuc_p = residue_two['O5\''].coord #phosphate position, DNA residues don't carry the ["P"] atom for some reason
        else:
            nuc_p = residue_two['P'].coord #phosphate position for RNA bases
   

        atom_name = nucleic_base_atoms.get(residue_two.get_resname()) #N1 or N9 for atom position
        nuc_b = residue_two[atom_name].coord #finding coordinate of that atom in the residue
        diff_vector_p = res_one_coord - nuc_p #subtracting vectors
        diff_vector_b = res_one_coord - nuc_b


    if is_amino_acid(residue_two): # if residue two is the amino acid

        res_two_coord = residue_two["CA"].coord #CA coordinate
        if residue_one.get_resname() in DNA_bases:
            if "P" in residue_two:
                nuc_p = residue_two['P'].coord
            else:
                nuc_p = residue_two['O5\''].coord #approximate phosphate position for DNA bases
        else:
            nuc_p = residue_one['P'].coord #phosphate position
        atom_name = nucleic_base_atoms.get(residue_two.get_resname())
        nuc_b = residue_one[atom_name].coord
        diff_vector_p = res_two_coord - nuc_p
        diff_vector_b = res_two_coord - nuc_b


    return np.sqrt(np.sum(diff_vector_p * diff_vector_p)), np.sqrt(np.sum(diff_vector_b * diff_vector_b)) #returns two distances



def calc_nuc_nuc_dist(residue_one, residue_two): #bit repetitive with thing above, make a seperate method with the overlap(?)
    """calculates the distances between pairs of the phosphate backbone and base of two nucleic acids

    Args:
        residue_one (residue class, amino acid): must be amino acid
        residue_two (residue class, amino acid): must be amino acid

    Returns:
        one_p_two_p : distance between phosphate of both nucleic acids
        one_p_two_b: distance between phosphate of first nucleic acid and base of second residue
        one_b_two_p: distance between base of first nucleic acid and phosphate of second
        one_b_two_b: distance between base of both residues
    """

    nucleic_base_atoms = {
        "A": "N9", "G": "N9", "C": "N1", "T": "N1", "U": "N1",
        "DA": "N9", "DG": "N9", "DC": "N1", "DT": "N1", "DU": "N1"
    }

    DNA_bases = {"DA", "DT", "DG", "DC", "DU"}


    if residue_one.get_resname() in DNA_bases:
        # print("new DNA molecule, name is ", residue_one.get_resname())
        # for atom in residue_one:
        #     print("Dna atom", atom)
        if "P" in residue_two:
            nuc_one_p = residue_two['P'].coord
        else:
            nuc_one_p = residue_two['O5\''].coord
    else:
        nuc_one_p = residue_one['P'].coord #phosphate position

    if residue_two.get_resname() in DNA_bases:
        if "P" in residue_two:
            nuc_two_p = residue_two['P'].coord
        else:
            nuc_two_p = residue_two['O5\''].coord
    else:
        nuc_two_p = residue_one['P'].coord #phosphate position

    diff_vector_one_p_two_p = nuc_one_p - nuc_two_p #phosphate of both nucleic acids
    diff_vector_one_p_two_b = nuc_one_p - residue_two[nucleic_base_atoms.get(residue_two.get_resname())].coord #phosphate of first nucleic acid and base of second residue
    diff_vector_one_b_two_p = residue_one[nucleic_base_atoms.get(residue_one.get_resname())].coord - nuc_two_p #base of first nucleic acid and phosphate of second
    diff_vector_one_b_two_b = residue_one[nucleic_base_atoms.get(residue_one.get_resname())].coord - residue_two[nucleic_base_atoms.get(residue_two.get_resname())].coord #base of both residues
    #find way to declutter :(
    
    return np.sqrt(np.sum(diff_vector_one_p_two_p * diff_vector_one_p_two_p)), np.sqrt(np.sum(diff_vector_one_p_two_b * diff_vector_one_p_two_b)), np.sqrt(np.sum(diff_vector_one_b_two_p * diff_vector_one_b_two_p)), np.sqrt(np.sum(diff_vector_one_b_two_b * diff_vector_one_b_two_b))






def run_padding(chain, position, pad_distance): # 3 on each side
    """returns a padded sequence of residues around a residue

    Args:
        chain (chain class object): chain A, B, C, etc
        position (int): position of "central" residue that the padding is built around on either side 
        pad_distance (int): default = 3

    Returns:
        string: sequence with padding

    
    Example:

    chain A position 5 (arginine) is passed through, with default pad distance of 3

    returned sequence: "KGGALPS" --> A is in the center, with the 3 residues around it
    """
    

    if (position - pad_distance  + 1< 1): #if position is at the start of the chain
        residues = [chain[res_id] for res_id in range(0, position + pad_distance)]
    elif (position + pad_distance >= len(chain)): #if position is at the end of the chain
        residues = [chain[res_id] for res_id in range(position - pad_distance, len(chain))]
    else:
        residues = [chain[res_id] for res_id in range(position - pad_distance , position + pad_distance + 1)] 


    # Convert residues to one-letter codes and concatenate them into a string
    if len(residues[0].get_resname()) == 3:
        sequence = ''.join([protein_letters_3to1[(residue.get_resname())] for residue in residues])
    elif len(residues[0].get_resname()) == 2: #checking for if it's DNA, where bases are defined as DA, DC, DG, etc. then, .get_resname can be used instead of having to convert text formats
        sequence = ''.join([(residue.get_resname())[1] for residue in residues]) #using only the second character, ie C, G of DC, DG
    else: #only one character, means it's an RNA
        sequence = ''.join([(residue.get_resname()) for residue in residues])
    # print(f"position {position - pad_distance} to {position + pad_distance + 1} is {sequence}")
    return sequence


def calc_dist_matrix(chain_one, chain_two, distance_threshold, pad_distance):
    """iterates between all pairs of residues of these two chains and returns various lists of the residue pairs that have distances below the distance threshold

    Args:
        chain_one (chain class object): first chain
        chain_two (chain class object): second chain
        distance_threshold (int): threshold for distance between two residues, default = 10
        pad_distance (int): distance for padding on either side, default = 3

    Returns:
        first residues (list, str): list of the first residues that had a distance with another residue that was less than the distance threshold, each residue is formatted as its number then amino acid, for example 25GLY
        second residues (list, str): list of the corresponding second residues that had a distance with another residue that was less than the distance threshold
        distances (list, int): list of the distances (already filtered for < distance threshold) between the residue pairs
        first_pads (list, str): padded residue sequences for  the first residue (residue one) in a pair
        second_pads (list, str): padded residue sequences for  the second residue (residue two) in a pair
        first_residues_num (list, str): only the numerical position of the same residues in first_residues
        second_residues_num (list, str): only the numerical position of the same residues in second_residues


    """
    first_residues = []
    second_residues = []
    distances = []
    first_pads = []
    second_pads = []
    first_residues_num = [] #these lists are only numbers, not amino acids
    second_residues_num = []
    answer = np.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one) : #iterating through all residue pairs between the two chains

        for col, residue_two in enumerate(chain_two) :

            if (is_amino_acid(residue_one) and is_amino_acid(residue_two)):

                answer[row, col] = calc_amino_dist(residue_one, residue_two)
                distance = calc_amino_dist(residue_one, residue_two)

                if (distance < distance_threshold): #checking if the paired residue distance is less than the paramterized threshold
                    #if it is, then add to list of residues
                    first_residues.append(f"{row}{residue_one.get_resname()}")
                    second_residues.append(f"{col}{residue_two.get_resname()}")
                    #run the padding on neighboring residues
                    first_pads.append(run_padding(chain_one, row, pad_distance))
                    second_pads.append(run_padding(chain_two, col, pad_distance))
                    #list to keep track of purely numerical indexes for qualifying residues
                    first_residues_num.append(f"{row}")
                    second_residues_num.append(f"{col}")
                    
                    distances.append(distance)
            
            if (is_nucleic_acid(residue_one) and is_nucleic_acid(residue_two)): #distances between two nucleic acids (all combinations of their bases and backbones)


                one_p_two_p, one_p_two_b, one_b_two_p , one_b_two_b = calc_nuc_nuc_dist(residue_one, residue_two)

                if one_p_two_p < distance_threshold:
                    
                    first_residues.append(f"{row}P")
                    second_residues.append(f"{col}P")


                    first_pads.append(run_padding(chain_one, row, pad_distance))
                    second_pads.append(run_padding(chain_two, col, pad_distance))

                    first_residues_num.append(f"{row}")
                    second_residues_num.append(f"{col}")
                    
                    distances.append(one_p_two_p)

                
                if one_p_two_b < distance_threshold:
                    
                    first_residues.append(f"{row}P")
                    second_residues.append(f"{col}B")


                    first_pads.append(run_padding(chain_one, row, pad_distance))
                    second_pads.append(run_padding(chain_two, col, pad_distance))

                    first_residues_num.append(f"{row}")
                    second_residues_num.append(f"{col}")
                    
                    distances.append(one_p_two_b)


                if one_b_two_p < distance_threshold:
                    
                    first_residues.append(f"{row}B")
                    second_residues.append(f"{col}P")

                    first_pads.append(run_padding(chain_one, row, pad_distance))
                    second_pads.append(run_padding(chain_two, col, pad_distance))

                    first_residues_num.append(f"{row}")
                    second_residues_num.append(f"{col}")
                    
                    distances.append(one_b_two_p)

                if one_b_two_b < distance_threshold:
                    
                    first_residues.append(f"{row}B")
                    second_residues.append(f"{col}B")

                    
                    first_pads.append(run_padding(chain_one, row, pad_distance))
                    second_pads.append(run_padding(chain_two, col, pad_distance))

                    first_residues_num.append(f"{row}")
                    second_residues_num.append(f"{col}")

                    distances.append(one_b_two_b)
                



            elif (is_nucleic_acid(residue_one) or is_nucleic_acid(residue_two)): #distances between one amino acid and one nucleic acid (CA of amino acid to backbone and base of the nucleic acid)

                one_p, one_b = calc_amino_nuc_dist(residue_one, residue_two) 

                
                if (is_nucleic_acid(residue_one)):

                    if one_p < distance_threshold:

                        first_residues.append(f"{row}P")
                        second_residues.append(f"{col}{residue_two.get_resname()}")

                        first_pads.append(run_padding(chain_one, row, pad_distance))
                        second_pads.append(run_padding(chain_two, col, pad_distance))

                        first_residues_num.append(f"{row }")
                        second_residues_num.append(f"{col}")
                        
                        distances.append(one_p)

                    if one_b < distance_threshold:

                        first_residues.append(f"{row}B")
                        second_residues.append(f"{col}{residue_two.get_resname()}")

                        first_pads.append(run_padding(chain_one, row, pad_distance))
                        second_pads.append(run_padding(chain_two, col, pad_distance))

                        first_residues_num.append(f"{row}")
                        second_residues_num.append(f"{col}")
                        
                        distances.append(one_p)
                
                if (is_nucleic_acid(residue_two)):

                    if one_p < distance_threshold:

                        first_residues.append(f"{row}{residue_one.get_resname()}")
                        second_residues.append(f"{col}P")

                        first_pads.append(run_padding(chain_one, row, pad_distance))
                        second_pads.append(run_padding(chain_two, col, pad_distance))

                        first_residues_num.append(f"{row }")
                        second_residues_num.append(f"{col}")
                        
                        distances.append(one_p)

                    if one_b < distance_threshold:


                        first_residues.append(f"{row}{residue_one.get_resname()}")
                        second_residues.append(f"{col}B")

                        first_pads.append(run_padding(chain_one, row, pad_distance))
                        second_pads.append(run_padding(chain_two, col, pad_distance))

                        first_residues_num.append(f"{row}")
                        second_residues_num.append(f"{col}")
                        
                        distances.append(one_b)




    return first_residues, second_residues, answer, distances, first_pads, second_pads, first_residues_num, second_residues_num






def create_final_columns(first_residues_num, first_pads):
    """synthesizes all of the padded ranges into final "binding interfaces"

    Args:
        first_residues_num (list): list of residue indicies 
        first_pads (list): list with all the padded sequences

    Returns:
        final_ranges (list): list containing "final strings" that represent the synthesized ranges

    note: the final ranges are synthesized by combining any of the padded residue ranges that are continuous. for example, if the residues 8, 9, and 10 
    have a residue pair < the distance threshold with sequences ANNL, NNLG, and NLGG, the final range would be ANNLGG

    """
    for i, res in enumerate(first_residues_num):
        first_residues_num[i] = int(res)

    final_ranges = []
    current_string = ""
    for i, res in enumerate(first_residues_num): 
        if i == 0:
            num_consec = 1 #keeping track of the # of consecutive residues
            current = first_residues_num[0]
            current_string = first_pads[0]


        if i >= len(first_residues_num) - 1: #if we have reached the end of the list of residues 

            final_ranges.append(current_string) #append whatever current string there is so far       
            return final_ranges


        if first_residues_num[i] == first_residues_num[i + 1]: #if the next residue is equal to the current residue (for example a residue might have <10 A with 2 different residues on the other chain)
            num_consec += 1
        elif first_residues_num[i] == first_residues_num[i+1] - 1: ##directly above
            if len(first_pads[i]) <= len(first_pads[i+1]):
                num_consec = 1
                current += 1
                current_string = current_string + (first_pads[i+1])[-1] 
            else:
                num_consec = 1
                current += 1
        else: ## big jump for example 1 --> 10, where it should be considered as seperate final ranges

            num_consec = 1 #resetting
            current = first_residues_num[i+1]
            final_ranges.append(current_string)
            # print("final string", current_string)
            current_string = first_pads[i+1]

def calculate_similarity(score, seq2):
    """given the "score" of the best alignment, returns the percentage of bases that are an exact match (what the score is score)


    Args:
        score (int): score of the best alignment generated by the aligner -> alignments[0]
        seq2 (string): shorter sequence that is being aligned to the binding site

    Returns:
        _type_: _description_
    """

    length = (len(seq2))
    return (score / length) * 100

def sequence_alignments(binding_sites, seq2): #seq1 is longer sequence, seq2 is shorter sequence
    """given the binding site and comparison sequence, 

    Args:
        binding_sites (_type_): _description_
        seq2 (_type_): _description_

    Returns:
        _type_: _description_
    """

    scores = []
    aligner = Align.PairwiseAligner() #creating the aligner from Biopython

    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    

    for seq1 in binding_sites:
        alignments = aligner.align(seq1, seq2)


    # Print the alignment in a nice format
    # for alignment in alignments:
    #     print(format_alignment(*alignment))


    # Calculate similarity for the best alignment
        best_alignment = alignments[0]
        # print("best alignment", best_alignment)

        similarity_percentage = calculate_similarity(best_alignment.score, seq2)
        scores.append(similarity_percentage)

        # print(f"Percentage Sequence Similarity: {similarity_percentage:.2f}%")

    return max(scores)


def identity_score(binding_sites, seq2):
    # determine which sequence is longer
    scores = []
    for seq1 in binding_sites:
        if len(seq1) < len(seq2):
            seq1, seq2 = seq2, seq1
        
        # print("seq 1 seq2", seq1, seq2)

        
        max_score = float('-inf')
        best_alignment = None

        # Slide the shorter sequence along the longer sequence
        for i in range(len(seq1) - len(seq2) + 1):
            score = 0
            for j in range(len(seq2)):
                if seq1[i + j] == seq2[j]:
                    score += 1  # Each match contributes a score of 1

            if score > max_score:
                max_score = score
                best_alignment = i  
        # return max_score
            
        scores.append(max_score)
            
        # print("max score", max_score)
    # print("scores", scores)
    if (max(scores) / len(seq2)) * 100 > 100:
        return 100
    else:
        return (max(scores) / len(seq2)) * 100

def aa_residues(chain):
    """removes residues that aren't a nucleic or amino acid

    Args:
        chain (chain class object): chain A, B, etc

    Returns:
        (list): list of only the amino acid + nucleic acid residues
    """

    nucleic_acids = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"}

    aa_only = []
    for i in chain:
        if (i.get_resname() in standard_aa_names) or (i.get_resname() in nucleic_acids):
            aa_only.append(i)
    return aa_only



#####################

def extract_structure(file_path, file_type):
    """given a file path and type, it uses the corresponding parser to extract the structure object

    Args:
        file_path (str): path to pdb/cif file
        file_type (_type_): _description_

    Raises:
        ValueError: file type besides cif or pdb

    Returns:
        structure class object: _description_
    """
    if file_type == 'cif':
        parser = MMCIFParser(QUIET=True)
    elif file_type == 'pdb':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file type: only 'cif' and 'pdb' are supported.")
    
    structure_id = os.path.splitext(os.path.basename(file_path))[0]
    structure = parser.get_structure(structure_id, file_path)
    return structure

def create_csv(model, chain_number_1, chain_number_2, pdb_name, distance_threshold, pad_distance, binding_sites):
    """ puts all of the elements together to create a csv between two chains

    Args:
        model (model class): subunit of a structure
        chain_number_1 (chain class): chain 1
        chain_number_2 (chain class): chain 2
        pdb_name (string): pdb name, for example "5tph"
        distance_threshold (int): threshold for residue distances
        pad_distance (int): distance for padding on both sides
    """

    d = dict(enumerate(string.ascii_uppercase, 1))
    ch1 = (d[chain_number_1]) #converts the number to string
    ch2 = (d[chain_number_2])


    chain_a = aa_residues(model[ch1])
    chain_b = aa_residues(model[ch2])

    first_residues, second_residues, dist_matrix, distances, first_pads, second_pads, first_residues_num, second_residues_num = calc_dist_matrix(chain_a, chain_b, distance_threshold, pad_distance)


    for i in range (len(second_residues_num)):
        second_residues_num[i] = int(second_residues_num[i])


    zipped_lists = list(zip(second_residues_num, second_pads))

    sorted_zipped_lists = sorted(zipped_lists, key = lambda x: x[0])

    res = [[i for i, j in sorted_zipped_lists],
       [j for i, j in sorted_zipped_lists]]

    second_residues_num = list(res[0])
    second_pads_ = list(res[1])

    total_length = len(first_pads)

    final_ranges_one = (create_final_columns(first_residues_num, first_pads))

    sim_scores_1 =  []
    id_scores_1 = []
    for bind_range in final_ranges_one:
        sim_scores_1.append(round(sequence_alignments(binding_sites[f"Sequences Chain {ch1}"], bind_range), 2))
        id_scores_1.append(round(identity_score(binding_sites[f"Sequences Chain {ch1}"], bind_range), 2))

    # print("sim scores", sim_scores)
    sim_scores_1_copy = sim_scores_1[:]
    id_scores_1_copy = id_scores_1[:]
    final_ranges_1_copy = final_ranges_one[:]


    if final_ranges_one != None:
        if len(final_ranges_one) < total_length:
            for i in range (total_length - len(final_ranges_one)):
                final_ranges_one.append("")
                sim_scores_1.append("")
                id_scores_1.append("")



    #final ranges for second residue in pairs

    final_ranges_two = (create_final_columns(second_residues_num, second_pads_))


    sim_scores_2 = []
    id_scores_2 = []

    for bind_range in final_ranges_two:
        sim_scores_2.append(round(sequence_alignments(binding_sites[f"Sequences Chain {ch2}"], bind_range), 2))
        id_scores_2.append(round(identity_score(binding_sites[f"Sequences Chain {ch2}"], bind_range), 2))

    sim_scores_2_copy = sim_scores_2[:]
    id_scores_2_copy = id_scores_2[:]
    final_ranges_2_copy = final_ranges_two[:]

    if final_ranges_two != None:
        if len(final_ranges_two) < total_length:
            for i in range (total_length - len(final_ranges_two)):
                final_ranges_two.append("") #pandas only accepts arrays of equal length in the dataframe, so add blank elements to reach equal length
                sim_scores_2.append("")
                id_scores_2.append("")
        
    chain_one = [f"{ch1}"] * len(first_residues)
    chain_two= [f"{ch2}"] * len(first_residues)



    # print("final string length", len(final_ranges_one))


    # print ("Minimum distance", numpy.min(dist_matrix))
    # print ("Maximum distance", numpy.max(dist_matrix))

    # import pylab
    # pylab.matshow(numpy.transpose(dist_matrix))
    # pylab.colorbar()
    # pylab.show()



    max_length = max(len(first_residues), len(second_residues))
    pdb_names = [f"{pdb_name}"] * max_length

    main_data = {
        "pdb name": pdb_names,
        "chain one": chain_one,
        "chain two": chain_two,
        "first residue": first_residues,
        "second residue": second_residues,
        "distance between": distances,
        "first residue padding": first_pads,
        "second residue padding": second_pads,
        "final ranges for residue one": final_ranges_one,
        "final ranges for residue two": final_ranges_two,
        "similarity scores for res 1": sim_scores_1,
        "similarity scores for res 2": sim_scores_2
    }


    data_scores = {
        "final ranges for residue one": final_ranges_one,
        "similarity scores for res 1": sim_scores_1,
        "identity scores res1": id_scores_1,
        "final ranges for residue two": final_ranges_two,
        "similarity scores for res 2": sim_scores_2,
        "identity scores res2": id_scores_2

    }



    df = pd.DataFrame(main_data)
    df_2 = pd.DataFrame(data_scores)



    print("data frame between chains", ch1, ch2)
    print (df)


    return df, df_2, sim_scores_1_copy , id_scores_1_copy, final_ranges_1_copy, sim_scores_2_copy , id_scores_2_copy, final_ranges_2_copy 




# root = tk.Tk()
# root.withdraw()
# data_file_path = filedialog.askopenfilename()
# #specifiying file paths and names
# output_file_path = (os.path.split(data_file_path))[0]
# data_file_name = (os.path.split(data_file_path))[1]
def main(input_dir, output_dir, binding_file_name, distance_threshold, pad_distance):
    """_summary_

    Args:
        input_dir (string): _description_
        output_dir (string): _description_
        distance_threshold (float): _description_
        pad_distance (int): _description_
    """

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)

    # Ensure input directory exists
    if not os.path.isdir(input_dir):
        print(f"Error: The input directory '{input_dir}' does not exist.")
        return
    
    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Find all .cif and .pdb files in the input directory
    cif_files = glob.glob(os.path.join(input_dir, '*.cif'))
    pdb_files = glob.glob(os.path.join(input_dir, '*.pdb'))
    
    all_files = cif_files + pdb_files

    columns = defaultdict(list) # each value in each column is appended to a list

    # os.chdir(input_dir)

    # print(os.getcwd())
    with open(binding_file_name) as f:
        reader = csv.DictReader(f) # read rows into a dictionary format
        for row in reader: # read a row as {column1: value1, column2: value2,...}
            for (k,v) in row.items(): # go over each column name and value 
                columns[k].append(v) # append the value into the appropriate list
    
    # os.chdir(output_dir)


    binding_sites = columns
    final_sequence_scores = []


    if not all_files:
        print(f"No .cif or .pdb files found in the directory '{input_dir}'.")
        return
    
    ranges_1 = []
    sim_scores_1= []
    id_scores_1 = []
    ranges_2 = []
    sim_scores_2 = []
    id_scores_2 = []
    file_names = []
    chain_names = []

    all_lists = []

    d = dict(enumerate(string.ascii_uppercase, 1))

    # Process each file
    for file_path in all_files:
        file_type = os.path.splitext(file_path)[1][1:]  # get the file extension without the dot
        structure = extract_structure(file_path, file_type)
        pdb_name = os.path.splitext(os.path.basename(file_path))[0]
        data_frames = []

        model = structure[0]
        num_chains = len(list(model))

        print("Processing file:", file_path)
        print("Number of chains:", num_chains)
        for i in range(num_chains):
            for j in range(i + 1, num_chains):
                current_df, df_2, sim_scores_1_copy , id_scores_1_copy, final_ranges_1_copy, sim_scores_2_copy , id_scores_2_copy, final_ranges_2_copy  = create_csv(structure[0], i + 1, j + 1, pdb_name, distance_threshold, pad_distance, columns)
                data_frames.append(current_df)
                ranges_1.extend(final_ranges_1_copy)
                ranges_1.extend(final_ranges_2_copy)
                sim_scores_1.extend(sim_scores_1_copy)
                sim_scores_1.extend(sim_scores_2_copy)
                id_scores_1.extend(id_scores_1_copy)
                id_scores_1.extend(id_scores_2_copy)
                file_names.extend([f"{pdb_name}"] * (len(final_ranges_1_copy) + len(final_ranges_2_copy)))
                chain_names.extend([d[i + 1]]*len(final_ranges_1_copy)) 
                chain_names.extend([d[j + 1]]*len(final_ranges_2_copy)) 



        
        result = pd.concat(data_frames)
        result.to_csv(f"{pdb_name}_binding_interface.csv",  index=False)



    all_lists.append(file_names)
    all_lists.append(chain_names)
    all_lists.append(ranges_1)
    all_lists.append(id_scores_1)
    all_lists.append(sim_scores_1)


    max_length = max(len(lst) for lst in all_lists)

    # Pad lists with empty strings
    padded_lists = [lst + [""] * (max_length - len(lst)) for lst in all_lists]

    headers = ["file_name", "chain", "residue 1 ranges", "res 1 id score", "res 1 sim score"] #, "res 2 ranges", "res 2 id score", "res 2 sim scores"

    sequence_data = pd.DataFrame({f'col_{i+1}': padded_lists[i] for i in range(len(padded_lists))})

    sequence_data.columns = headers

    sequence_data = sequence_data.sort_values(by="res 1 id score", ascending=False)




    sequence_data.to_csv(f"final_sequences.csv",  index=False)






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process .pdb files in a directory.')
    
    parser.add_argument('--input_dir', type=str, required=True, help='The input file directory containing .pdb files')
    parser.add_argument('--output_dir', type=str, required=True, help='The output file directory')
    parser.add_argument('--binding_file_name', type = str, required = True, help = 'CSV with list of binding sites (sequences)')
    parser.add_argument('--distance_threshold', type=float, default = 10, required=False, help='distance threshold in angstroms')
    parser.add_argument('--padding', type=int, default = 3, required=False, help='# of residues added as "padding" on either side of the residues')

    args = parser.parse_args()
    
    main(args.input_dir, args.output_dir, args.binding_file_name, args.distance_threshold, args.padding)



