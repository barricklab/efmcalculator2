#!/usr/bin/env python

import unittest
import os
import subprocess
import csv
import shutil
import pathlib
import efmcalculator
import Bio.SeqIO as SeqIO


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
output_directory_path = os.path.join(THIS_DIR)


def csv_are_identical(csv_file_path_1, csv_file_path_2):
    'Utility function for comparing output, Returns whether files are identical as bool'
    failed = False
    try:
        with open(csv_file_path_1) as csv1:
            reader = csv.reader(csv1)
            csv_1_contents = [row for row in reader]
    except:
        print(f'Failed to open {str(csv_file_path_1)}')
        failed = True

    try:
        with open(csv_file_path_2) as csv2:
            reader = csv.reader(csv2)
            csv_2_contents = [row for row in reader]
    except:
        print(f'Failed to open {str(csv_file_path_2)}')
        failed = True

    if len(csv_1_contents) != len(csv_2_contents):
        print(f'CSVs are different lengths')
        failed = True

    differences = []
    for i in range(0, len(csv_1_contents)):
        if csv_1_contents[i] != csv_2_contents[i]:
            # print(csv_1_rows[i], "\n!=\n",csv_1_rows[i])
            differences.append(i)
    if differences:
        print(f'CSVs have {len(differences)} different values at rows {differences}')
        failed = True

    if failed == True:
        return False
    else:
        return True

##############################################################################################
# Unit tests (function calls)
##############################################################################################
class test_unit_run_complete_wo_errors(unittest.TestCase):
    """Tests to ensure that EFM2 is capable of completing a single example. Fails on error only"""
    def test_unit_run_complete_wo_errors(self):
        L6_10_plasmid = output_directory_path  + "/../examples/L6-10_plasmid_bba.gb"
        L6_10_plasmid = pathlib.Path(L6_10_plasmid)
        L6_10_plasmid = efmcalculator.parse_file(L6_10_plasmid)
        inseqs = list()
        inseqs.append(L6_10_plasmid)
        results = efmcalculator.predict_many(L6_10_plasmid, "pairwise", True)
        for result, seq in zip(results, inseqs):
            seq = seq[0]
            post_result = efmcalculator.post_process(*result, len(seq.seq), True)
if __name__ == "__main__":
    unittest.main()
