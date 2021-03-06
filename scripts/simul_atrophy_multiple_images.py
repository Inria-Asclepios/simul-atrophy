#!/usr/bin/env python
""" Script to run simul@atrophy for multiple images
"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
# Asclepios INRIA Sophia Antipolis
#

import os
import os.path as op
import argparse as ag
import bish_utils as bu
import time as tm

def get_input_options():
    """ command line interface, get input options and interact with the user
    """
    #parse the command line inputs
    parser = ag.ArgumentParser(
        'Read list of patient ids from a file and for of those patients launch '
        'job in cluster to run AdLemModel with the options given in ops_file.')
    parser.add_argument(
        'simdir_list', help='File with list of directories present at '
        '$ADLEM_DIR/results/patients/ to simulate images for. Each of these '
        ' directories must have all the inputs specified in the -ops_file.')
    parser.add_argument(
        'ops_file', help='File containing all the options required by '
        'runAdLemModel.py except the patients_list.')
    parser.add_argument('--print_only', action='store_true', help='just print'
                        ', dont do anything!')
    ops = parser.parse_args()
    return ops

def main():
    """ Runs the script. Run with -h to see options.
    """
    ops = get_input_options()
    adlem_dir = os.getenv('ADLEM_DIR')
    run_adlem = 'python ' + op.join(adlem_dir, 'scripts/runAdLemModel.py')
    if adlem_dir is None:
        raise ValueError('environment variable ADLEM_DIR not set.')
    adlem_ops = bu.get_lines_as_list(ops.ops_file, '#')
    sim_ids = bu.get_lines_as_list(ops.simdir_list, '#')
    for sim_id in sim_ids:
        cmd = ('%s %s %s\n'
               % (run_adlem, sim_id, ' '.join(adlem_ops)))
        if ops.print_only:
            print cmd
        else:
            bu.print_and_execute(cmd)
            tm.sleep(0.1)

if __name__ == "__main__":
    main()
