#!/home/methyl/anaconda3/bin/python
"""
output_summary

Used to write output summary file that will be accessible to the user
"""
import json
import pandas
import argparse

import os
import csv
from collections import defaultdict

def read_signature_labels():

    label_map = {}
    if os.path.exists('signature_label.txt'):
        signaturelabelfilehandler = open('signature_label.txt', 'r')
        for line in signaturelabelfilehandler:
            line = line.rstrip()
            fields = line.split(' ')
            label_map[fields[1]] = fields[0]
        signaturelabelfilehandler.close()
    return label_map

def get_actual_label(label_map, key):

    if key in label_map:
        return label_map[key]
    return ''

def consolidate_geminals(options):

    opt_dict = {}
    new_options = []
    for o in options:
        resid = o.split('.')[0]
        if resid not in opt_dict:
            opt_dict[resid] = 1
            new_options.append(resid)

    return new_options

def activate_stereospecificity():

    stereospecific = False
    stereofilehandler = open("stereo.txt", "r")
    for line in stereofilehandler:
        line = line.rstrip()
        stereospecific = line
    stereofilehandler.close()

    if stereospecific == 'True':
        return "yes"
    else:
        return "no"

def generate_histogram(assignment_options):

    counter= {}
    counter[1] = 0
    counter[2] = 0
    counter[3] = 0
    counter[4] = 0
    counter[5] = 0
    histogram = []
    for a in assignment_options:
        options = a.split(' ')
        new_options = consolidate_geminals(options)
        if len(new_options) == 1:
            counter[1] += 1
        elif len(new_options) == 2:
            counter[2] += 1
        elif len(new_options) == 3:
            counter[3] += 1
        elif len(new_options) == 4:
            counter[4] += 1
        elif len(new_options) > 4 :
            counter[5] += 1


    for key in counter:
        k = '='+str(key)
        if key == 5:
            k = '>4'
        histogram.append(f"[no. of valid assignment options{k}]:{counter[key]:<4}| "+ counter[key] * "\u25a7")

    return histogram

def get_input_option(key):

    input_file = "input.json"
    data = {}
    if os.path.exists(input_file):
        with open(input_file) as json_file:
            data = json.load(json_file)
            if key in data:
                return data[key]

    return None


def generate_output_assignment_labels(signature, labels, sc_label, label_map):

    options = labels.split(' ')
    new_options = ''
    opts_dict = {}
    if labels:
        for label in options:
            text = list(label)
            color = text[0]
            resid = label.split('.')
            new_label = resid[0]
            atom = ''
            if color == 'V':
                if "?" not in signature:
                    if resid[1] == '1':
                        atom = 'R'
                    elif resid[1] == '2':
                        atom = 'S'
                else:
                    atom = '?'
                    if sc_label == 'pror':
                        atom = 'R'
                    elif sc_label == 'pros':
                        atom = 'S'
            elif color == 'L':
                if "?" not in signature:
                    if resid[1] == '1':
                        atom = 'R'
                    elif resid[1] == '2':
                        atom = 'S'
                else:
                    atom = '?'
                    if sc_label == 'pror':
                        atom = 'R'
                    elif sc_label == 'pros':
                        atom = 'S'

            opt = new_label+atom
            if opt not in opts_dict:
                opts_dict[opt] = 1
                if new_options:
                    new_options += ' '+opt
                else:
                    new_options += opt


    return (new_options)

def read_resonances():

    label_map = read_signature_labels()

    sc_label = get_input_option('sclabeling')

    resonances = []
    assignment_options = []
    if os.path.exists("assigned_2d.csv"):
        with open("assigned_2d.csv") as resonancefilehandler:
            csv_reader = csv.reader(resonancefilehandler, delimiter=',')
            title = True
            no_of_rows = 0
            for row in csv_reader:
                if not title:
                    if no_of_rows == 6:
                        r0 = get_actual_label(label_map, row[0])
                        asgs = row[2].split(' ')
                        r2 = get_actual_label(label_map, asgs[0])
                        for i in range(1, len(asgs)):
                            a_label = get_actual_label(label_map, asgs[i])
                            r2 += ' '+a_label
                        r5 = generate_output_assignment_labels(r0, row[5], sc_label, label_map)
                        assignment_options.append(row[5])
                        resonances.append('{:20}{:15}{:30}{:10}{:10}{:30}'.format(r0, row[1], r2, row[3], row[4], r5))
                    elif no_of_rows == 5:
                        r0 = get_actual_label(label_map, row[0])
                        r4 = generate_output_assignment_labels(r0, row[4], sc_label, label_map)
                        assignment_options.append(row[4])
                        resonances.append('{:20}{:15}{:10}{:10}{:30}'.format(r0, row[1], row[2], row[3], r4))
                else:
                    no_of_rows = len(row)
                    if len(row) == 6:
                        resonances.append('{:20}{:15}{:30}{:10}{:10}{:30}\n'.format('Label', 'Residue Type', 'Forced assignments', 'Carbon', 'Hydrogen', 'Valid assignment options'))
                    elif len(row) == 5:
                        resonances.append('{:20}{:15}{:10}{:10}{:30}\n'.format('Label', 'Residue Type', 'Carbon', 'Hydrogen', 'Valid assignment options'))
                    title = False

            resonancefilehandler.close()

    return (resonances, assignment_options)


def get_total_unused_noes(unused_list):

    short = 0
    long = 0
    for u in unused_list:
        if 'short' in u:
            short += 1
        else:
            long += 1

    return (short, long)

def check_alignment_file():

    alignmentfile = 'align.txt'
    found_error = False
    align_string = ''
    if os.path.exists(alignmentfile):
        alignmentfilehandler = open(alignmentfile, 'r')
        for line in alignmentfilehandler:
            if found_error:
                align_string += line
            if 'error' in line:
                found_error = True
                if 'pdb' in line:
                    align_string = "Input PDB file has missing atoms in methyl bearing residues. Please fix them and rerun.\n"
                else:
                    align_string = "The sequences used in the NMR construct and the input PDB structure have mismatches. Please see the alignment strings below to find out mismatches.\n\n"
        alignmentfilehandler.close()

    return (found_error, align_string)

def read_geminal_noes():

    geminalnoehandler = open('geminal_noes.txt')
    geminal = defaultdict(list)
    for line in geminalnoehandler:
        line = line.rstrip()
        fields = line.split(' ')
        if fields[0] in geminal:
            if fields[1] not in geminal[fields[0]]:
                geminal[fields[0]].append(fields[1])
        else:
            geminal[fields[0]].append(fields[1])

        if fields[1] in geminal:
            if fields[0] not in geminal[fields[1]]:
                geminal[fields[1]].append(fields[0])
        else:
            geminal[fields[1]].append(fields[0])

    geminalnoehandler.close()
    return geminal

def noe_class_count(noe_list, geminal_noes):

    no_symm = 0
    no_symm_gem = 0
    no_complex = 0
    simple_2 = 0
    simple_3 = 0
    for u in noe_list:
        if int(noe_list[u]) == 1:
            if u not in geminal_noes:
                no_symm += 1
            else:
                no_symm_gem += 1
        elif int(noe_list[u]) == 2:
            simple_2 += 1
        elif int(noe_list[u]) == 3:
            simple_3 += 1
        elif int(noe_list[u]) > 3:
                no_complex += 1

    return (no_symm, no_symm_gem, no_complex, simple_2, simple_3)

def print_crosspeak_information(outputfilehandler, diagonal_list, crosspeaks, noe_before_list, noe_after_list, geminal_noes):

    outputfilehandler.write("NOE peaks information\n========================\n")
    outputfilehandler.write("Diagonal NOEs\n===============\n")
    for i in range(0, len(diagonal_list)):
        outputfilehandler.write(diagonal_list[i]+"\n")

    outputfilehandler.write("\nNon-diagonal NOEs\n===============\n")
    for i in range(0, len(crosspeaks)):
        crosspeaks[i] = crosspeaks[i].strip();
        crosspeaks[i] = crosspeaks[i].replace("\s+", " ")
        fields = crosspeaks[i].split(" ")
        label = fields[0]
        if fields[0] in noe_after_list:
            noe_class = ''
            if int(noe_after_list[fields[0]]) == 1:
                if fields[0] in geminal_noes:
                    noe_class = "geminal cluster"
                else:
                    if int(noe_before_list[fields[0]]) > 1:
                        noe_class = "No symmetry matches after reduction"
                    else:
                        noe_class = "No symmetry matches"
            elif int(noe_after_list[fields[0]]) > 3:
                noe_class = "Complex - not used"
            outputfilehandler.write(crosspeaks[i]+"    "+noe_class+"\n")
        else:
            outputfilehandler.write(crosspeaks[i]+"\n")

    outputfilehandler.write("---------------------------------------------------------------------------------------\n\n")


def main():

    (align_error, align_string) = check_alignment_file()
    if not align_error:

        short_noes = 0
        long_noes = 0
        short_noes_diagonal = 0
        long_noes_diagonal = 0
        used_noes_total_long = 0
        used_noes_total_short = 0
        used_noes_total_diagonal = 0
        peaks_2d = 0
        used_short = 0
        used_long = 0
        used_short_aft_clustering = 0
        used_long_aft_clustering = 0
        exclusive = 0
        isolated = 0

        crosspeaks = []
        crosspeak_found = False

        geminals = []
        geminals_found = False

        ground_truth = []
        ground = False

        noe_after_list = {}
        noe_after_tag = False

        noe_before_list = {}
        noe_before_tag = False

        diagonal_list = []
        diagonal_found = False

        histogram = []

        dc = ''
        edc = ''

        annotationfilehandler = open("noe_annotated.txt", "r")
        for line in annotationfilehandler:
            line = line.rstrip()
            if "end after list" in line:
                noe_after_tag = False
            if "end before list" in line:
                noe_before_tag = False
            if noe_after_tag == True:
                fields = line.split(" ")
                noe_after_list[fields[0]] = fields[1]
            if noe_before_tag == True:
                fields = line.split(" ")
                noe_before_list[fields[0]] = fields[1]
            if "end crosspeaks" in line:
                crosspeak_found = False
            if crosspeak_found == True:
                crosspeaks.append(line)
            if "end geminals" in line:
                geminals_found = False
            if "end ground truth" in line:
                ground = False
            if "end diagonals" in line:
                diagonal_found = False
            if geminals_found == True:
                geminals.append(line)
            if ground == True:
                ground_truth.append(line)
            if diagonal_found == True:
                diagonal_list.append(line)
            if "effective no. of assignment options" in line:
                histogram.append(line)
            if "processed" in line:
                if "peaks" in line:
                    used_noes_total = line.split(":")[1]
                elif "diagonal" in line:
                    used_noes_total_diagonal = line.split(":")[1]
            if "long.csv" in line:
                if "peaks" in line:
                    long_noes = line.split(":")[1]
                elif "diagonal" in line:
                    long_noes_diagonal = line.split(":")[1]
            if "short.csv" in line:
                if "peaks" in line:
                    short_noes = line.split(":")[1]
                elif "diagonal" in line:
                    short_noes_diagonal = line.split(":")[1]
            if "2D" in line:
                peaks_2d = line.split(":")[1]
            if "used" in line:
                if "short" in line:
                    if "clustering" in line:
                        used_short_aft_clustering = line.split(":")[1]
                    else:
                        used_short = line.split(":")[1]
                if "long" in line:
                    if "clustering" in line:
                        used_long_aft_clustering = line.split(":")[1]
                    else:
                        used_long = line.split(":")[1]
            if "Exclusively" in line:
                exclusive = line.split(":")[1]
            if "Isolated" in line:
                isolated = line.split(":")[1]
            if "start crosspeaks" in line:
                crosspeak_found = True
            if "Geminal pairs" in line:
                geminals_found = True
            if "start ground truth" in line:
                ground = True
            if "start after list" in line:
                noe_after_tag = True
            if "start before list" in line:
                noe_before_tag = True
            if "start diagonals" in line:
                diagonal_found = True
            if "degree_connectivity" in line and "effective_degree_connectivity" not in line:
                dc = line.split(" ")[1]
            if "effective_degree_connectivity" in line:
                edc = line.split(" ")[1]


        outputfilehandler = open("final_output.txt", "w")
        outputfilehandler.write(f"Thank you for using MAUS. If you have any questions or face any issues with this software, please e-mail methyl@ucsc.edu.\n\n")

        outputfilehandler.write(f"Below is the input data you submitted:\n")
        outputfilehandler.write(f"PDB file:{get_input_option('pdbfile')}\n")
        outputfilehandler.write(f"PDB chain:{get_input_option('pdbchain')}\n")
        outputfilehandler.write(f"PDB file type:{get_input_option('pdb_file_type')}\n")
        outputfilehandler.write(f"Labeling scheme:{get_input_option('labeling')}\n")
        outputfilehandler.write(f"Side chain labeling scheme:{get_input_option('sclabeling')}\n")

        if get_input_option('pdb_file_type') == 'crystal':
            outputfilehandler.write(f"Relax option:{get_input_option('relax')}\n")

        outputfilehandler.write(f"2D HMQC reference file:{get_input_option('referencefile')}\n")
        if 'manual_rt' in get_input_option('residue_type'):
            outputfilehandler.write(f"Residue type/s:manual\n")
        else:
            outputfilehandler.write(f"Residue type/s:auto\n")
        if 'manual_g' in get_input_option('geminals'):
            outputfilehandler.write(f"Geminal NOEs:manual\n")
        else:
            outputfilehandler.write(f"Geminal NOEs:auto\n")
        outputfilehandler.write(f"Known assignments:{get_input_option('known_assignments')}\n")

        if 'dimension-3d' in get_input_option('dimension'):
            outputfilehandler.write(f"Long mixing time NOE file:{get_input_option('longnoe_3d')}\n")
            outputfilehandler.write(f"Short mixing time NOE file:{get_input_option('shortnoe_3d')}\n")
            outputfilehandler.write(f"Maximum NOE distance for long-range NOEs:{get_input_option('radius_3d')}\n")
            outputfilehandler.write(f"Maximum NOE distance for short-range NOEs:{get_input_option('shortradius_3d')}\n")
            outputfilehandler.write(f"Diagonal NOEs tolerances:{get_input_option('diag_tol_C_3d')}\n")
            outputfilehandler.write(f"Symmetry tolerances:{get_input_option('sym_tol_C_3d')}\n")
            outputfilehandler.write(f"Clustering tolerances:{get_input_option('clustering_tol_C_3d')}, {get_input_option('clustering_tol_H_3d')}\n")
        else:
            outputfilehandler.write(f"NOE file:{get_input_option('longnoe_4d')}\n")
            outputfilehandler.write(f"Maximum NOE distance:{get_input_option('radius_4d')}\n")
            outputfilehandler.write(f"Diagonal NOEs tolerances:{get_input_option('diag_tol_C_4d')}, {get_input_option('diag_tol_H_4d')}\n")
            outputfilehandler.write(f"Symmetry tolerances:{get_input_option('sym_tol_C_4d')}, {get_input_option('sym_tol_H_4d')}\n")
            outputfilehandler.write(f"Clustering tolerances:{get_input_option('clustering_tol_C_4d')}, {get_input_option('clustering_tol_H_4d')}\n")

        outputfilehandler.write(f"Are the input assignments classified as stereospecific?:{activate_stereospecificity()}\n")

        outputfilehandler.write("\n")

        (resonances, assignment_options) = read_resonances()
        geminal_noes = read_geminal_noes()
        if len(resonances):

            (no_symm_bef, no_symm_gem_bef, no_complex_bef, simple_2_bef, simple_3_bef) = noe_class_count(noe_before_list, geminal_noes)
            (no_symm_aft, no_symm_gem_aft, no_complex_aft, simple_2_aft, simple_3_aft) = noe_class_count(noe_after_list, geminal_noes)

            total_aft = simple_2_aft+simple_3_aft

            outputfilehandler.write(f"Data-to-constraints summary\n======================\n")
            outputfilehandler.write(f"# 2D peaks: {peaks_2d}\n")
            outputfilehandler.write(f"Total NOEs: {len(crosspeaks) - 1 + len(diagonal_list)}\n")
            outputfilehandler.write(f"# diagonal NOEs: {len(diagonal_list)}\n")
            outputfilehandler.write(f"# NOEs without symmetry matches before/after reduction of complex components: {no_symm_bef}/{no_symm_aft}\n")
            outputfilehandler.write(f"# NOEs in geminal clusters: {no_symm_gem_aft}\n")
            outputfilehandler.write(f"# NOEs from simple components of size 2 before/after reduction of complex components: {simple_2_bef}/{simple_2_aft}\n")
            outputfilehandler.write(f"# NOEs from simple components of size 3 before/after reduction of complex components: {simple_3_bef}/{simple_3_aft}\n")
            outputfilehandler.write(f"# NOEs that form complex components before/after reduction of complex components: {no_complex_bef}/{no_complex_aft}\n")
            outputfilehandler.write(f"# NOEs used by MAUS: {total_aft}\n")
            #outputfilehandler.write(f"Degree connectivity: {dc}\n")
            outputfilehandler.write(f"Effective degree connectivity: {edc}\n\n")

            # NOEs used before and after reduction of complex components
            # left over NOEs

            print_crosspeak_information(outputfilehandler, diagonal_list, crosspeaks, noe_before_list, noe_after_list, geminal_noes)

            outputfilehandler.write("Resonance assignments statistics (at the residue level)\n===========================================================\n")

            # overwriting old histogram here.
            histogram = generate_histogram(assignment_options)

            for h in histogram:
                outputfilehandler.write(h+"\n")

            outputfilehandler.write("\n")
            outputfilehandler.write("Resonance assignments\n==============================\n")
            for i in range(0, len(resonances)):
                outputfilehandler.write(resonances[i]+"\n")

            outputfilehandler.write("---------------------------------------------------------------------------------------\n\n")

        else:
            outputfilehandler.write("MAUS did not successfully complete your job. Below are a few ways to troubleshoot and rerun MAUS.\n\n")
            outputfilehandler.write("1. In order to ensure that NOEs are clustered to their correct 2D methyl resonances, the 3D or 4D spectra must be carefully aligned, i.e., the 3D or 4D NOE data have to be phase corrected. In addition, the NOE peaks must be picked at sufficiently high signal-to-noise ratio levels (>5).\n")
            outputfilehandler.write("2. If you have manually specified peak residue types, please check them to make sure you are 100% confident or use our residue type classifier by selecting \"Auto\" under residue type classifier in the input form.\n")
            outputfilehandler.write("3. If you have manually specified geminal methyl peaks, please check them to make sure you are 100% confident or use our geminal classifier by selecting \"Auto\" under geminal classifier in the input form.\n")
            outputfilehandler.write("4. If you have provided known assignments and forced them, please check to make sure you are 100% confident about the assignments.\n")
            outputfilehandler.write("5. Try increasing the clustering tolerances in the input form. If your clustering tolerances are very low, then an NOE peak may not be clustered into the correct 2D peak/s.\n")
            outputfilehandler.write("6. Try increasing the symmetry tolerances in the input form. If your symmetry tolerances are very low, then your NOE network may not contain the correct symmetry connectivities.\n\n")

            outputfilehandler.write("Below is the NOE network (or H) contructed using the information provided. Please check if the peaks cluster and symmetrize correctly. The NOE peaks that do not have any symmetry matches or have ambiguous symmetry matches are not used in the construction of H.\n\n")
            print_crosspeak_information(outputfilehandler, diagonal_list, crosspeaks, noe_before_list, noe_after_list, geminal_noes)

            if len(ground_truth) > 0:
                outputfilehandler.write("Please find below an additional message that shows why your NOE network could not fit into input structure within the specified radii:\n\n")
                for i in range(0, len(ground_truth)):
                    outputfilehandler.write(ground_truth[i]+"\n")


        outputfilehandler.close()
        annotationfilehandler.close()

    else:
        outputfilehandler = open("final_output.txt", "w")
        outputfilehandler.write(f"Thank you for using MAUS. If you have any questions or face any issues with this software, please e-mail methyl@ucsc.edu.\n\n")
        outputfilehandler.write(align_string)
        outputfilehandler.close()


if __name__ == "__main__":

    main()
