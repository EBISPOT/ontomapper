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
# import spotilities as spot
from spotilities import newsflash
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
    9. oxo_url       : URL of OxO web service; default = public OxO server
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
        newsflash("Spreadsheet location is %s" % spreadsheet)
        newsflash("Column index is %d" % colno)

        url_bool = re.compile('[a-zA-Z](\w|[-+.])*://.*')
        # filestuff = None
        if url_bool.match(spreadsheet):
            t0 = time.time()
            newsflash("Getting spreadsheet from URL ...")
            r = requests.get(spreadsheet, allow_redirects=True)
            t1 = time.time()
            newsflash("It took %f seconds to retrieve the spreadsheet." % float(t1 - t0))
            filestuff = io.StringIO(r.content.decode('utf-8'))
        else:
            filestuff = spreadsheet

        newsflash("Pandafying spreadsheet ...")
        source_df = pd.read_table(filestuff, low_memory=False, keep_default_na=False)
        # source_df = source_df.head(6000)

        newsflash("Getting source terms ...")
        source_iris = source_df.iloc[:, colno]
        # iri_lists = source_iris.apply(lambda x: x.split(", "))
        """
        1. Return empty array if empty string.
        2. Might want to split by comma only, then remove spaces: more robust!
        """
        newsflash("Breaking source term strings into lists ...")
        iri_lists = source_iris.apply(lambda x: [] if x == '' else x.split(", "))
        """ Use nested lambdas to collect source IRIs from Series of lists """
        newsflash("Generate dictionary keyed on unique source terms ...")
        iri_lists.apply(lambda x, y: list(map(lambda z: y.update({z: None}), x)), args=[iri_map])
        newsflash("Dictionary generated!")
        ss_dict.update({'pandafued': source_df})
        ss_dict.update({'unique_iris': iri_map})
        newsflash("Returning from function 'parse_ss' ...")
    except requests.exceptions.InvalidSchema as is_error:
        """ pandas should have coped with distinguishing text file from URL already """
        # if is_error.
        #     raise
        newsflash('Error retrieving spreadsheet?')
        newsflash(is_error)
        # raise
    return ss_dict


def map_iris(iri_dict, target_ontologies, threshold, use_paxo, oxo_inner_url, query_size, verbose):
    quantified_url = "%s?size=%d" % (oxo_inner_url, query_size)
    data = {'ids': list(iri_dict.keys()), 'mappingTarget': target_ontologies}
    json_strings = []
    data['distance'] = '1' if threshold >= 100 else '2' if threshold >= 75 else '3' if threshold >= 50 else ''
    """ If boundary less than 50%, throw 'confidence too low' error: need to code! """
    oxo_hit_counter = 0
    while quantified_url is not None:
        reply = requests.post(quantified_url, data)
        oxo_hit_counter += 1
        json_content = reply.content
        json_string = json.loads(json_content)
        newsflash(json_string, verbose)
        json_strings.append(json_string)
        try:
            these_results = json_string["_embedded"]["searchResults"]
        except:
            newsflash("IRI map load aborted: OxO hit %d times, with %d query terms" %
                           (oxo_hit_counter, oxo_hit_counter * query_size))
            raise
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
            # newsflash(ontology_dict)
        try:
            quantified_url = json_string["_links"]["next"]["href"]
        except KeyError:
            newsflash("Stopped")
            quantified_url = None
    newsflash("No. of iterative calls to OxO web API was %d" % len(json_strings))
    """ Passed-in dictionary object is mutated in situ: no need to return it """
    return None


def augment(panda_input, iri_map, table_format, colno, keep_original, iri_format):
    colname = panda_input.columns[colno]
    prev_colname = panda_input.columns[colno - 1]
    next_colname = panda_input.columns[colno + 1]
    out_columns = panda_input.columns
    newsflash('Here are the pandas columns')
    newsflash(out_columns)
    newsflash()
    in_tuple_counter = 0
    # out_tuple_counter = 0
    out_dict_list = []
    if table_format in {'uni-column', 'multi-column'}:
        extra_col_dict_list = []
    tt0 = time.time()
    tt1 = tt0
    for in_tuple in panda_input.itertuples():
        """ Need to convert back to regular tuple, from pandafied named tuple with extra leading index number """
        in_supple = tuple(in_tuple[1:])
        # source_string = in_tuple[colname]
        source_string = in_supple[colno]
        source_terms = source_string.split(", ")
        target_groups = {}
        if table_format == 'in-situ' and keep_original:
            """ Key '00source00' is lazy, collational way of placing source terms at top of list, where we want them """
            target_groups.setdefault('00source00', source_terms )
        for source_term in source_terms:
            # newsflash("Source term is %s" % source_term)
            # map_dict = iri_map[source_term]
            map_dict = iri_map.get(source_term)
            if map_dict:
                for m in map_dict:
                    """ target_groups assigned one key per target ontology --- NOT per source term in source cell! """
                    target_groups.setdefault(m, []).append(map_dict[m])
                    # newsflash("Found term %s" % map_dict[m])
        for target_group in target_groups:
            target_groups[target_group] = ', '.join(target_groups[target_group])
        tg_series = pd.Series(target_groups)
        # tg_series = tg_series.reindex(sorted(tg_series.index))
        # newsflash(tg_series)

        # out_dict_list = []

        if table_format in {'in-situ', 'uni-row', 'uni-column'}:
            target_string = ', '.join(tg_series.values)
            # newsflash("Additions to %s from MeSH: %s" % (source_string, target_string))
            # newsflash()

        if table_format in {'in-situ', 'uni-row', 'multi-row'}:
            out_dict = dict(zip(out_columns, in_supple))
            # newsflash(out_dict)

            if table_format == 'in-situ' or keep_original:
                if len(tg_series) > 0 or keep_original:
                    if table_format == 'in-situ':
                        out_dict[colname] = target_string
                    out_dict_list.append(out_dict)

            if table_format == 'uni-row' and len(tg_series) > 0:
                out_dict_extra = dict(out_dict)
                out_dict_extra[colname] = target_string
                out_dict_list.append(out_dict_extra)

            elif table_format == 'multi-row':
                for hit in tg_series.values:
                    out_dict_iter = dict(out_dict)
                    out_dict_iter[colname] = hit
                    out_dict_list.append(out_dict_iter)

        elif table_format in {'uni-column', 'multi-column'} and (len(tg_series) > 0 or keep_original):
            """ Now need to handle uni- and multi-column outputs """
            out_dict = dict(zip(out_columns, in_supple))
            out_dict_list.append(out_dict)
            extra_col_dict = {}
            if table_format == 'multi-column':
                extra_col_dict = dict(tg_series)
            elif table_format == 'uni-column':
                # extra_col_dict = dict({'all_ontologies', ', '.join(tg_series.values)})
                extra_col_dict = dict({colname: target_string})
            extra_col_dict_list.append(extra_col_dict)

        in_tuple_counter += 1
        if in_tuple_counter % 4000 == 0:
            tt2 = time.time()
            newsflash("Processed %d thousand input records: took %.2f s (increment of %.2f s)" %
                           (int(in_tuple_counter / 1000), float(tt2 - tt0), float(tt2 - tt1)))
            tt1 = tt2

    newsflash("No. of records in output spreadsheet is %d" % len(out_dict_list))
    panda_output = pd.DataFrame(out_dict_list, columns=out_columns)
    if table_format in {'uni-column', 'multi-column'}:
        newsflash("Adding new columns ...")
        extra_columns_df = pd.DataFrame(extra_col_dict_list,
                                        columns=[colname] if table_format == 'uni-column' else tg_series.keys())
        panda_output = pd.concat([panda_output.loc[:, :colname if keep_original else prev_colname], extra_columns_df,
                                  panda_output.loc[:, next_colname:]], axis=1)

    return panda_output


def re_ontologise(input, output, format, column_index, keep, target, uri_format, boundary, paxo, oxo_url, quantity, verbose):

    target = sorted(target)
    newsflash("Length of target ontology array is %d" % len(target))
    for t in target:
        newsflash("Target is %s" % t)
    ss_dict = parse_ss(input, column_index)
    iri_map = ss_dict['unique_iris']
    """ Print out list of source iris in iri_map """
    # iri_counter = 0
    # for src_iri in iri_map:
        ### print("%d\t%s\t%s" % (iri_counter, src_iri, iri_map[src_iri]))  # Print values _and_ keys
        # newsflash("%d\t%s" % (iri_counter, src_iri))
        # iri_counter += 1
    panda_original = ss_dict['pandafued']
    newsflash("Calling map_iris with url = '%s' ..." % oxo_url)
    map_iris(iri_map, target, boundary, paxo, oxo_url, quantity, verbose)
    newsflash("Calling augment ...")
    gwas_enriched = augment(panda_original, iri_map, format, column_index, keep, uri_format)
    """ Print out augmented_panda here ... """
    newsflash("No. of dictionary elements: %d" % len(ss_dict))
    newsflash("No. of rows in spreadsheet: %d" % len(panda_original))
    newsflash("No. of unique IRIs: %d" % len(iri_map))
    newsflash('', verbose)
    newsflash(ss_dict['unique_iris'], verbose)
    # newsflash(ss_dict.keys())
    # for spot_key in ss_dict:
    #     newsflash(spot_key)
    """ Enable print to check for uniqueness of IRIs """
    # for iri_key in ss_dict['unique_iris']:
    #     newsflash(iri_key)
    # newsflash(ss_dict['unique_iris'])
    newsflash("Outputting enriched GWAS spreadsheet ...")
    # print(gwas_enriched.head(30).to_csv(index=False, sep='\t'))
    print(gwas_enriched.to_csv(index=False, sep='\t'))


def main():

    """ First of all, check configuration file """
    parser = argparse.ArgumentParser(description='Type of URI to output, etc.')
    parser.add_argument('-g', '--config', default='./ontomapper.ini')

    """ Second of all, get configuration info from configuration file """
    ontoconfig = configparser.ConfigParser(os.environ)
    ontoconfig.optionxform = str
    namespace, extra = parser.parse_known_args()
    ontoconfig.read(namespace.config)

    """ Third of all, parse the rest of the switches, possibly using defaults from configuration file """
    parser.add_argument('-i', '--input', default=ontoconfig.get('Params', 'gwas_spreadsheet'))
    parser.add_argument('-o', '--output', default='')
    parser.add_argument('-f', '--format', default='multi-column',
                        choices=['in-situ', 'uni-column', 'multi-column', 'uni-row', 'multi-row'])
    parser.add_argument('-c', '--column-index', type=int, default=35)
    kmeg = parser.add_mutually_exclusive_group(required=False)
    kmeg.add_argument('-k', '--keep', dest='keep', action='store_true')
    kmeg.add_argument('-d', '--no-keep', dest='keep', action='store_false')
    parser.set_defaults(keep=True)
    parser.add_argument('-t', '--target', default=['mesh'], nargs='+')
    parser.add_argument('-u', '--uri-format', default='long')
    parser.add_argument('-b', '--boundary', type=int, default=100)
    parser.add_argument('-x', '--oxo_url', default=ontoconfig.get('Params', 'oxo_url'))
    parser.add_argument('-p', '--paxo', type=bool, nargs='?', const=True, default=False)
    parser.add_argument('-q', '--quantity', type=int, default=ontoconfig.getint('Params', 'api_record_quantity'))
    parser.add_argument('-v', '--verbose', type=bool, nargs='?', const=True, default=False)
    args = parser.parse_args()

    """ vars returns a dictionary from the Namespace object; """
    arg_dict = vars(args)
    arg_dict.pop('config')
    newsflash(arg_dict, arg_dict['verbose'])

    """ ** unpacks a dictionary """
    re_ontologise(**arg_dict)


if __name__ == "__main__":
    main()
