#!/usr/bin/env python3
"""Check that all the input files were processed and produced an output with
the expected number of events for each input file."""

import os
import argparse

import h5py
import numpy as np
from ROOT import TFile # pylint: disable=E0611
from larcv import larcv # pylint: disable=W0611


def main(source, source_list, output, dest, suffix):
    """Checks the output of the SPINE process.

    The script loops over the input files, check that there is an output file
    in the expected location and further checks that the output file entry
    count matches that of the input file.

    Produces a list of input files that have no or incomplete output in a text
    file (the name of which is provided with the `-o` or `--output` flag. This
    can be used to reprocess missing/incomplete input files.

    .. code-block:: bash

        $ python3 -c SPINE_CONFIG -S missing_list.txt

    Parameters
    ----------
    source : List[str]
        List of paths to the input files
    source_list : str
        Path to a text file containing a list of data file paths
    output : str
        Path to the output text file with the list of badly processed files
    dest : str
        Destination directory for the original SPINE process
    suffix : str
        Suffix added to the end of the input files by the original SPINE process
    """
    # If using source list, read it in
    if source_list is not None:
        with open(source_list, 'r', encoding='utf-8') as f:
            source = f.read().splitlines()

    # Initialize the output text file
    out_file = open(output, 'w', encoding='utf-8')

    # Loop over the list of files in the input
    print("\nChecking existence and completeness of output files.")
    miss_list, inc_list = [], []
    for idx, file_path in enumerate(source):
        # Find the base name of the input file (without extension)
        base = os.path.basename(file_path)
        stem, _ = os.path.splitext(base)

        # Check that the output exists under the expected path
        out_base = f'{stem}_{suffix}.h5'
        out_path = f'{dest}/{out_base}'
        if not os.path.isfile(out_path):
            print(f"- Missing: {out_base}")
            out_file.write(f'{file_path}\n')
            miss_list.append(file_path)
            break

        # If the output does exist, check that the input and output have the
        # same number of entries. First count the number of entries in the
        # first tree of the input file (assumed to all match, as they should)
        f = TFile(file_path)
        key = [key.GetName() for key in f.GetListOfKeys()][0]
        num_entries = f.Get(key).GetEntries()
        f.Close()

        # Then check the number of events in the output file
        with h5py.File(out_path) as f:
            if len(f['events']) != num_entries:
                print(f"-Incomplete: {out_base}")
                out_file.write(f'{file_path}\n')
                inc_list.append(file_path)

    num_miss = len(miss_list)
    num_inc = len(inc_list)
    print(f"\nFound {num_miss + num_inc} problematic output file(s):")
    print(f"- {num_miss} missing output file(s);")
    print(f"- {num_inc} incomplete output file(s).")

    # Close text file
    out_file.close()


if __name__ == "__main__":
    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description="Check dataset validity")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--source', '-s',
                       help='Path or list of paths to data files',
                       type=str, nargs="+")
    group.add_argument('--source-list', '-S',
                       help='Path to a text file of data file paths',
                       type=str)

    parser.add_argument('--output', '-o',
                        help='Path to the output file',
                        type=str)

    parser.add_argument('--dest',
                        help='Destination directory for the original SPINE process',
                        type=str)
        
    parser.add_argument('--suffix',
                        help='Suffix added to the input files by the original SPINE process',
                        type=str)

    args = parser.parse_args()

    # Execute the main function
    main(args.source, args.source_list, args.output, args.dest, args.suffix)
