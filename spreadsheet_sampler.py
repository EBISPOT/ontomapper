#!/usr/bin/env python3

import argparse
import io
import configparser
import os
import pandas as pd
import re
import requests
from spotilities import newsflash
from spotilities import config_or_bust
import sys
import time


"""
__PARAMETERS__ (to main program: all non-positional key/value pairs, to be initialised with default values):

    1. input-file    : URL or filename; default = online GWAS spreadsheet
    2. output        : filename; default = standard output
    3. file-format   : input / output spreadsheet file format: 'csv' or 'tsv'; default = 'tsv'
    4. sample-size   : count of records required in output sample
"""


def sample_ss(input_file, output, file_format, sample_size, column_name):
    ss_dict = {}
    iri_map = {}
    try:
        # newsflash("Spreadsheet location is %s" % spreadsheet)
        # newsflash("Column index is %d" % colno)
        separator = "\t" if file_format == 'tsv' else ","

        url_bool = re.compile('[a-zA-Z](\w|[-+.])*://.*')
        # filestuff = None
        if url_bool.match(input_file):
            t0 = time.time()
            newsflash("Getting spreadsheet from URL ...")
            r = requests.get(input_file, allow_redirects=True)
            t1 = time.time()
            newsflash("It took %.2f seconds to retrieve the spreadsheet" % float(t1 - t0))
            filestuff = io.StringIO(r.content.decode('utf-8'))
        else:
            filestuff = input_file

        newsflash("Pandafying spreadsheet ...")
        source_df = pd.read_csv(filestuff, sep=separator, low_memory=False, keep_default_na=False)

        newsflash("Generating random sample of records ...")
        output_df = source_df.sample(n=sample_size).loc[:, column_name]
        print(output_df.to_csv(index=False, sep=separator))

    except requests.exceptions.InvalidSchema as is_error:
        """ pandas should have coped with distinguishing text file from URL already """
        newsflash('Error retrieving spreadsheet?')
        newsflash(is_error)
    return None


def main():

    """ First of all, check configuration file """
    parser1 = argparse.ArgumentParser(description='Config filepath', add_help=False)
    parser1.add_argument('-g', '--config', help='filepath of config file (ini format)')
    namespace, extra = parser1.parse_known_args()
    config_file = namespace.config

    """ Second of all, get configuration info from configuration file (if specified) """
    ontoconfig = configparser.ConfigParser(os.environ)
    ontoconfig.optionxform = str
    column_list = None
    if config_file is not None:
        ontoconfig.read(config_file)
        columns_plus = ontoconfig.items('Columns')
        spurious_columns = ontoconfig.items('DEFAULT')
        real_columns = dict(set(columns_plus) - set(spurious_columns))
        column_list = list(real_columns.values())

    """ Third of all, parse the rest of the switches, possibly using defaults from configuration file """
    parser2 = argparse.ArgumentParser(prog='Ontomapper',
                                      description="%s%s" %
                                                   ('Allows the generation of speadsheet subsample. Specifically, ',
                                                    'samples rows at random and columns by parameterised request.'),
                                      parents=[parser1])
    cfg_sect_lookup = config_or_bust(ontoconfig, 'Params')
    parser2.add_argument('-i', '--input-file', default=cfg_sect_lookup('input_file', 'string'),
                         help='location of input spreadsheet: accepts filepath or URL')
    parser2.add_argument('-o', '--output', default=cfg_sect_lookup('output', 'string'),
                         help='output spreadsheet filepath **NO CURRENT EFFECT**')
    parser2.add_argument('-f', '--file-format', choices=['csv', 'tsv'],
                         default=cfg_sect_lookup('file_format', 'string'),
                         help='file format (both input and output)')
    parser2.add_argument('-s', '--sample-size', type=int, default=cfg_sect_lookup('sample_size', 'int'),
                         help='number of records to return in randomly sampled spreadsheet')
    parser2.add_argument('-c', '--column-name', nargs='+', default=column_list,
                         help='space-separated list of column names required in sampled output')

    if len(sys.argv) < 2:
        parser2.print_help(sys.stderr)
        newsflash()
        sys.exit(0)

    args = parser2.parse_args()

    """ vars returns a dictionary from the Namespace object; """
    arg_dict = vars(args)
    newsflash()
    newsflash("These are your opening parameters:")
    newsflash()
    # newsflash("Length of args dictionary is %d" % len(arg_dict))
    newsflash(pd.Series(arg_dict))
    newsflash()

    """ config is still in arg_dict at this point """
    arg_dict.pop('config')

    """ Don't check values of reserved options, which have no effect at the moment """
    active_arg_dict = arg_dict.copy()
    for inactive_arg in ['output']:
        active_arg_dict.pop(inactive_arg)
    if None in active_arg_dict.values():
        newsflash()
        newsflash("Please set values for the following parameters---on command line or in config file!")
        newsflash()
        for cfg_key in active_arg_dict:
            if active_arg_dict[cfg_key] is None:
                newsflash("\t%s" % cfg_key)
        newsflash()
        sys.exit(1)

    """ '**' unpacks a dictionary """
    sample_ss(**arg_dict)


if __name__ == "__main__":
    main()
