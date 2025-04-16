#!/usr/bin/env python

import unittest
import os
import subprocess
import csv
import shutil
import pathlib
import efmcalculator
import Bio.SeqIO as SeqIO
import polars as pl
import efmcalculator
from itertools import chain

from efmcalculator.ingest.bad_state_mitigation import detect_special_cases

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
        print("Test: Run complete without errors")
        L6_10_plasmid = output_directory_path  + "/../examples/1_L6-10_plasmid.gb"
        L6_10_plasmid = pathlib.Path(L6_10_plasmid)
        L6_10_plasmid = efmcalculator.ingest.parse_file(L6_10_plasmid)
        inseqs = list()
        statemachine = efmcalculator.StateMachine.StateMachine()
        statemachine.import_sequences(inseqs)
        for seqobject in statemachine.user_sequences.values():
            seqobject.call_predictions()


class test_unit_filtering_ssr(unittest.TestCase):
    """Tests SSR filtering. Fails when an SSR doesn't get filtered appropriately"""
    def test_unit_filtering_ssr(self):
        return
        test_sequence = "CACACACACA"+"ATGTTTTCACACACACA"*6 +"ATGTTTTGGGAGAGAGCACACAATGTTTT"
        shouldnt_exist = ["ACACACACAATGTTTTC",
                          "CACACACAATGTTTTCA",
                          "ACACACAATGTTTTCAC",
                          "CACACAATGTTTTCACA",
                          "ACACAATGTTTTCACAC",
                          "CACAATGTTTTCACACA",
                          "ACAATGTTTTCACACAC",
                          "CAATGTTTTCACACACA",
                          "AATGTTTTCACACACAC",
                          "ATGTTTTCACACACACA",
                          "AC",
                          "A",
                          "C",
                          "G",
                          "AG"]
        should_exist = ["CACACACACAATGTTTT"]

        for is_circular in [True, False]:
            if not is_circular:
                should_exist += ["T", "CA"]

            ssr_df, _, _ = efmcalculator.predict(test_sequence, "pairwise", is_circular)
            ssr_df = efmcalculator.pipeline.filtering.filter_ssrs(ssr_df, len(test_sequence), is_circular)

            ssr_df = ssr_df.with_columns(
                pl.lit(shouldnt_exist).alias("test_case")
            ).with_columns(
                pl.col("repeat").is_in(pl.col("test_case")).alias("fail")
            )
            #assert (True,) not in ssr_df.select("fail").unique().rows(), "SSR failing to filter some nested sequences"
            for positive_check in should_exist:
                detected_ssrs = ssr_df.filter(pl.col("repeat") == positive_check)
                assert detected_ssrs.height >= 1, "Known SSR is missing"

class test_unit_filtering_rmd(unittest.TestCase):
    def test_unit_filtering_rmd(self):
        test_sequence = "AGTCCAA" + "AAAAAAAAAAAAAAAAAAAAA" + "AGTCCAA"
        is_circular = False
        ssr_df, srs_df, rmd_df = efmcalculator.predict(test_sequence, "pairwise", is_circular)

if __name__ == "__main__":
    unittest.main()
