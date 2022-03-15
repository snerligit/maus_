#!/usr/bin/python


#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: October 18, 2021
#   Email: snerli@ucsc.edu
#

# import other required libraries
import os
import sys
import json
import argparse

params = {}
params["email"] = "null"
params["pdbfile"] = "null"
params["pdbchain"] = "null"
params["pdb_file_type"] = "nmr"
params["labeling"] = "ailv"
params["sclabeling"] = "none"
params["relax"] = "quick"
params["is_symmetric"] = "symmetric_no"
params["sequence"] = "null"
params["cross"] = "null"
params["referencefile"] = "null"
params["residue_type"] = "manual_rt"
params["geminals"] = "manual_g"
params["dimension"] = "dimension-3d"
params["radius_3d"] = "10"
params["shortradius_3d"] = "8"
params["longnoe_3d"] = "null"
params["shortnoe_3d"] = "null"
params["diag_tol_C_3d"] = "0.1"
params["sym_tol_C_3d"] = "0.15"
params["clustering_tol_C_3d"] = "0.15"
params["clustering_tol_H_3d"] = "0.02"
params["radius_4d"] =  "null"
params["longnoe_4d"] = "null"
params["diag_tol_C_4d"] = "null"
params["diag_tol_H_4d"] = "null"
params["sym_tol_C_4d"] = "null"
params["sym_tol_H_4d"] = "null"
params["clustering_tol_C_4d"] = "null"
params["clustering_tol_H_4d"] = "null"

def get_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Method to convert input text file to json file")
    parser.add_argument("-help", action='store_true', help='print input options')
    parser.add_argument("-gen", action='store_true', help='generate sample file')
    parser.add_argument("-infile", help="Input text file")
    parser.add_argument("-outfile", help="Output file that is used to execute main", default="default.sh")
    parser.add_argument("-workingdir", help="Working directory path of the example", default='./')
    parser.add_argument("-exe", help="full path to the main executable", default="./Linux/main")

    args = parser.parse_args()

    return args

def print_help():

    print ("Please provide a text file with following fields")
    for key in params:
        print(key)

    print ("You can use \"python preprocess_input.py -gen\" to generate a sample file named default.txt. You can then edit this file and add clustom details")

def generate_inputfile():

    inputfilehandler = open("default.txt", "w")
    for key, value in params.items():
        inputfilehandler.write(key+" "+value+"\n")
    inputfilehandler.close()

def read_inputfile(args):

    inputfilehandler = open(args.infile, "r")
    for line in inputfilehandler:
        line = line.rstrip()
        fields = line.split()
        if fields[0] not in params:
            print ("Invalid keyword: ", fields[0], " Run with help option to see params supported")
        params[fields[0]] = fields[1]
    inputfilehandler.close()

def convert2sh(args):

    #/Users/snerli/Work/UCSC/Research/maus/code_rewrite/maus_v2/src/main /Users/snerli/Work/UCSC/Research/maus/code_rewrite/maus_v2/test_files/datasets/b2m A ailv 0 HB2M_2D.list false false 3 10.0 6.5 HB2M_300ms_3D.list HB2M_50ms_3D.list 0.1 0.01 0.15 0.02 0.15 0.02 3 false false false false false false false false HB2M_ensemble_01.pdb HB2M_ensemble_05.pdb HB2M_ensemble_09.pdb HB2M_ensemble_13.pdb HB2M_ensemble_17.pdb HB2M_ensemble_02.pdb HB2M_ensemble_06.pdb HB2M_ensemble_10.pdb HB2M_ensemble_14.pdb HB2M_ensemble_18.pdb HB2M_ensemble_03.pdb HB2M_ensemble_07.pdb HB2M_ensemble_11.pdb HB2M_ensemble_15.pdb HB2M_ensemble_19.pdb HB2M_ensemble_04.pdb HB2M_ensemble_08.pdb HB2M_ensemble_12.pdb HB2M_ensemble_16.pdb HB2M_ensemble_20.pdb
    #~/anaconda3/bin/python3 /Users/snerli/Work/UCSC/Research/maus/code_rewrite/maus_v2/scripts/output_summary.py

    inputfilehandler = open(args.outfile, "w")
    inputfilehandler.write(f"{args.exe} {args.workingdir} ")
    inputfilehandler.write(f"{params['pdbchain']} {params['labeling']}")
    inputfilehandler.close()


if __name__ == "__main__":

    args = get_args()
    if args.help:
        print_help()
    elif args.gen:
        generate_inputfile()
    else:
        read_inputfile(args)
        convert2sh(args)
