#!/usr/bin/env python3

import argparse
# import csv
import io
import json
# from configparser import ConfigParser
import configparser
import os
import pandas as pd
import re
import requests
import spotilities as spot
import time


"""
__PARAMETERS__ (to main program: all non-positional key/value pairs, to be initialised with default values):

    1. input         : URL or filename; default = online GWAS spreadsheet
    2. output        : filename; default = standard output
    3. format        : format of new entries in output spreadsheet: 'in-situ' (requires/implies --keep), 'uni-column',
                       'multi-column', 'uni-row' or 'multi-row'; default = 'multi-column'
    4. column-index  : column index (starting at 0) to parse for source IRIs; default based on header search?
    5. keep          : whether to retain column containing source IRIs in output; boolean, default = no
    6. target        : ontolog[y|ies] to search for equivalent terms; default = MeSH
    7. uri-format    : return format of target IRIs (curies, long form, etc.); default = long form
    8. boundary      : boundary value (%age) for confidence level; default = 100 (OxO distance 1)
    9. oxo           : 'pub' or 'dev'; default = 'pub'
    10. paxo         : use 'experimental' PAXO code?; boolean, default = no
    11. config       : configuration file path
    12. quantity     : no. of query terms to include in single API call
    13. verbose      : do we want a bunch of large json objects dumped to standard error?; boolean, default = no
    14. ????         : source ontology from which terms in existing spreadsheet come; default = 'EFO' --- Probably
                       unnecessary! --- check!!
    15. ????         : input format of source IRIs?
    16. ????         : negation operators for the three boolean parameters? ... probably unnecessary.
    17. ????         : column header as mutually exlusive alternative to column index? --- does this have implications
                       for processing of the column, inside and outside pandas?
    18. ????         : option to change format of incoming IRIs before re-output (implies --keep)


__FUNCTIONS__

1. Return a dictionary containing two objects:
     a) the contents of the original spreadsheet itself, as some kind of pandas thing?
     b) a list of unique URIs
   params: see parameters list above
   return: two-element dictionary, as detailed above
    
2. Map each URI in the original list to its equivalent (via OxO) in one or more other ontologies, returning a
     dictionary with original URI as key, and a nested dictionary keyed on ontology name, with equivalent URI
     in that ontology, plus either distance or normalised threshold incorporating distance, as the value.
   params: unique URI list
   return: dictionary

3. Return a modified spreadsheet based on the input speadsheet and the URI mappings
   params: a) pandas representation of original spreadsheet
           b) dictionary object containing URI mappings
   return: pandas representation of modified spreadsheet, including new ontology column(s)
"""


def parse_ss(spreadsheet, colno):
    ss_dict = {}
    iri_map = {}
    try:
        spot.newsflash("Spreadsheet location is %s" % spreadsheet)
        spot.newsflash("Column index is %d" % colno)

        url_bool = re.compile('[a-zA-Z](\w|[-+.])*://.*')
        # filestuff = None
        if url_bool.match(spreadsheet):
            t0 = time.time()
            spot.newsflash("Getting spreadsheet from URL ...")
            r = requests.get(spreadsheet, allow_redirects=True)
            t1 = time.time()
            spot.newsflash("It took %f seconds to retrieve the spreadsheet." % float(t1 - t0))
            filestuff = io.StringIO(r.content.decode('utf-8'))
        else:
            filestuff = spreadsheet

        spot.newsflash("Pandafying spreadsheet ...")
        source_df = pd.read_csv(filestuff, sep='\t', low_memory=False, keep_default_na=False)

        spot.newsflash("Getting source terms ...")
        source_iris = source_df.iloc[:, colno]
        # iri_lists = source_iris.apply(lambda x: x.split(", "))
        """
        1. Return empty array if empty string.
        2. Might want to split by comma only, then remove spaces: more robust!
        """
        spot.newsflash("Breaking source term strings into lists ...")
        iri_lists = source_iris.apply(lambda x: [] if x == '' else x.split(", "))
        """ Use nested lambdas to collect source IRIs from Series of lists """
        spot.newsflash("Generate dictionary keyed on unique source terms ...")
        iri_lists.apply(lambda x, y: list(map(lambda z: y.update({z: None}), x)), args=[iri_map])
        spot.newsflash("Dictionary generated!")
        ss_dict.update({'pandafued': source_df})
        ss_dict.update({'unique_iris': iri_map})
        spot.newsflash("Returning from function 'parse_ss' ...")
    except requests.exceptions.InvalidSchema as is_error:
        """ pandas should have coped with distinguishing text file from URL already """
        # if is_error.
        #     raise
        spot.newsflash('Error retrieving spreadsheet?')
        spot.newsflash(is_error)
        # raise
    return ss_dict


def map_iris(iri_dict, target_ontologies, threshold, use_paxo, oxo_inner_url, query_size, verbose):
    quantified_url = "%s?size=%d" % (oxo_inner_url, query_size)
    data = {'ids': list(iri_dict.keys()), 'mappingTarget': target_ontologies}
    json_strings = []
    data['distance'] = '1' if threshold >= 100 else '2' if threshold >= 75 else '3' if threshold >= 50 else ''
    """ If boundary less than 50%, throw 'confidence too low' error: need to code! """
    while quantified_url is not None:
        reply = requests.post(quantified_url, data)
        json_content = reply.content
        json_string = json.loads(json_content)
        spot.newsflash(json_string, verbose)
        json_strings.append(json_string)
        these_results = json_string["_embedded"]["searchResults"]
        for this_result in these_results:
            hits = this_result["mappingResponseList"]
            ontology_dict = {}
            for hit in hits:
                target_ontology = hit['targetPrefix']
                """ Create one key per target ontology, then append individual hits to associated array """
                ontology_dict.setdefault(target_ontology, []).append(hit['curie'])
            for okey in ontology_dict:
                ontology_dict[okey] = ', '.join(ontology_dict[okey])
            iri_dict[this_result["queryId"]] = ontology_dict
            # spot.newsflash(ontology_dict)
        try:
            quantified_url = json_string["_links"]["next"]["href"]
        except KeyError:
            spot.newsflash("Stopped")
            quantified_url = None
    spot.newsflash("No. of iterative calls to OxO web API was %d" % len(json_strings))
    """ Passed-in dictionary object is mutated in situ: no need to return it """
    return None


def augment(panda_input, iri_map, table_format, colno, keep_original, iri_format):
    colname = panda_input.columns[colno]
    out_columns = panda_input.columns
    spot.newsflash('Here are the pandas columns')
    spot.newsflash(out_columns)
    spot.newsflash()
    spot.newsflash('Here is the pandas dataframe')
    panda_output = pd.DataFrame(columns=out_columns)
    spot.newsflash(panda_output)
    for in_tuple in panda_input.itertuples():
        # source_string = in_tuple[colname]
        source_string = in_tuple[colno]
        source_terms = source_string.split(", ")
        target_groups = {}
        if table_format == 'in-situ':
            target_groups.setdefault('__source__', source_terms )
        for source_term in source_terms:
            # map_dict = iri_map[source_term]
            map_dict = iri_map.get(source_term)
            if map_dict:
                for m in map_dict:
                    """ target_groups assigned one key per target ontology --- NOT per source term in source cell! """
                    target_groups.setdefault(m, []).append(map_dict[m])
        for target_group in target_groups:
            target_groups[target_group] = ', '.join(target_groups[target_group])
        tg_series = pd.Series(target_groups)

        out_dict_list = []

        if table_format in {'in-situ', 'uni-row', 'uni-column'}:
            target_string = ', '.join(tg_series.values)

        if table_format in {'in-situ', 'uni-row', 'multi-row'}:
            out_dict = dict(zip(out_columns, in_tuple))

            if table_format == 'in-situ' or keep_original:
                if table_format == 'in-situ':
                    out_dict[colname] = target_string
                out_dict_list.append(out_dict)

            if table_format == 'uni-row':
                out_dict_extra = dict(out_dict)
                out_dict_extra[colname] = target_string
                out_dict_list.append(out_dict_extra)

            elif table_format == 'multi-row':
                for target, hit in tg_series:
                    out_dict_iter = dict(out_dict)
                    out_dict_iter[colname] = hit
                    out_dict_list.append(out_dict_iter)

        """ Now need to handle uni- and multi-column outputs """

        panda_output = panda_output.append(out_dict_list)

    return panda_output


def re_ontologise(input, output, format, column_index, keep, target, uri_format, boundary, paxo, oxo_url, quantity, verbose):

    spot.newsflash("Length of target ontology array is %d" % len(target))
    ss_dict = parse_ss(input, column_index)
    iri_map = ss_dict['unique_iris']
    panda_original = ss_dict['pandafued']
    spot.newsflash("Calling map_iris with url = '%s' ..." % oxo_url)
    map_iris(iri_map, target, boundary, paxo, oxo_url, quantity, verbose)
    spot.newsflash("Calling augment ...")
    gwas_enriched = augment(panda_original, iri_map, format, column_index, keep, uri_format)
    """ Print out augmented_panda here ... """
    spot.newsflash("No. of dictionary elements: %d" % len(ss_dict))
    spot.newsflash("No. of rows in spreadsheet: %d" % len(iri_map))
    spot.newsflash("No. of unique IRIs: %d" % len(panda_original))
    spot.newsflash('', verbose)
    spot.newsflash(ss_dict['unique_iris'], verbose)
    # spot.newsflash(ss_dict.keys())
    # for spot_key in ss_dict:
    #     spot.newsflash(spot_key)
    """ Enable print to check for uniqueness of IRIs """
    # for iri_key in ss_dict['unique_iris']:
    #     spot.newsflash(iri_key)
    # spot.newsflash(ss_dict['unique_iris'])
    spot.newsflash("Outputting enriched GWAS spreadsheet ...")
    print(gwas_enriched.to_csv())


def main():

    parser = argparse.ArgumentParser(description='Type of URI to output, etc.')
    parser.add_argument('-g', '--config', default='./ontomapper.ini')

    """ First of all, get configuration info from file """
    ontoconfig = configparser.ConfigParser(os.environ)
    ontoconfig.optionxform = str
    namespace, extra = parser.parse_known_args()
    ontoconfig.read(namespace.config)

    # parser.add_argument('-i', '--input', default=ontoconfig.get('Params', 'gwas_file'))
    parser.add_argument('-i', '--input', default=ontoconfig.get('Params', 'gwas_spreadsheet'))
    parser.add_argument('-o', '--output', default='')
    parser.add_argument('-f', '--format', default='multi-column',
                        choices=['in-situ', 'uni-column', 'multi-column', 'uni-row', 'multi-row'])
    parser.add_argument('-c', '--column-index', type=int, default=36)
    parser.add_argument('-k', '--keep', type=bool, nargs='?', const=True, default=False)
    parser.add_argument('-t', '--target', default=['mesh'], nargs='+')
    parser.add_argument('-u', '--uri-format', default='long')
    parser.add_argument('-b', '--boundary', type=int, default=100)
    parser.add_argument('-x', '--oxo', default='pub', choices=['pub', 'dev'])
    parser.add_argument('-p', '--paxo', type=bool, nargs='?', const=True, default=False)
    parser.add_argument('-q', '--quantity', type=int, default=ontoconfig.getint('Params', 'api_record_quantity'))
    parser.add_argument('-v', '--verbose', type=bool, nargs='?', const=True, default=False)
    args = parser.parse_args()
    """ vars returns a dictionary from the Namespace object; """
    arg_dict = vars(args)
    arg_dict.pop('config')
    spot.newsflash(arg_dict, arg_dict['verbose'])

    oxo_url = ''
    try:
        oxo_url = ontoconfig.get('Params', 'oxo_' + arg_dict['oxo'])
    except:
        spot.newsflash('OxO config type has no corresponding URL')
    arg_dict.pop('oxo')
    arg_dict.update({'oxo_url': oxo_url})

    spot.newsflash(arg_dict, arg_dict['verbose'])

    """ ** unpacks a dictionary """
    re_ontologise(**arg_dict)


if __name__ == "__main__":
    main()
