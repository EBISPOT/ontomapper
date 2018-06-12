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
from spotilities import config_or_bust
import sys
import time


"""
__PARAMETERS__ (to main program: all non-positional key/value pairs, to be initialised with default values):

    1. input-file    : URL or filename; default = online GWAS spreadsheet
    2. output        : filename; default = standard output
    3. layout        : layout of new entries in output spreadsheet: 'in-situ', 'uni-column', 'multi-column', 'uni-row'
                       or 'multi-row'; default = 'multi-column'
    4. file-format   : input / output spreadsheet file format: 'csv' or 'tsv'; default = 'tsv'
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


def parse_ss(spreadsheet, separator, colno):
    ss_dict = {}
    iri_map = {}
    try:
        # newsflash("Spreadsheet location is %s" % spreadsheet)
        # newsflash("Column index is %d" % colno)

        url_bool = re.compile('[a-zA-Z](\w|[-+.])*://.*')
        # filestuff = None
        if url_bool.match(spreadsheet):
            t0 = time.time()
            newsflash("Getting spreadsheet from URL ...")
            r = requests.get(spreadsheet, allow_redirects=True)
            t1 = time.time()
            newsflash("It took %.2f seconds to retrieve the spreadsheet" % float(t1 - t0))
            filestuff = io.StringIO(r.content.decode('utf-8'))
        else:
            filestuff = spreadsheet

        newsflash("Pandafying spreadsheet ...")
        source_df = pd.read_csv(filestuff, sep=separator, low_memory=False, keep_default_na=False)
        # source_df = source_df.head(6000)

        newsflash("Getting source terms ...")
        source_iris = source_df.iloc[:, colno]
        # iri_lists = source_iris.apply(lambda x: x.split(", "))
        """
        1. Return empty array if empty string.
        2. Split by comma only, then remove spaces afterwards: more robust!
        """
        newsflash("Breaking source term strings into lists ...")
        # iri_lists = source_iris.apply(lambda x: [] if x == '' else x.split(", "))
        iri_lists = source_iris.apply(lambda x: [] if x == '' else list(map(lambda w: w.strip(), x.split(","))))
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
            newsflash("Stopped: all good!")
            quantified_url = None
    newsflash("No. of iterative calls to OxO web API was %d" % len(json_strings))
    """ Passed-in dictionary object is mutated in situ: no need to return it """
    return None


def augment(panda_input, iri_map, table_layout, colno, keep_original, iri_format):
    colname = panda_input.columns[colno]
    prev_colname = panda_input.columns[colno - 1]
    next_colname = panda_input.columns[colno + 1]
    out_columns = panda_input.columns
    newsflash()
    newsflash('These are the column headers of your pandas DataFrame:')
    newsflash()
    newsflash(out_columns)
    newsflash()
    in_tuple_counter = 0
    # out_tuple_counter = 0
    out_dict_list = []
    if table_layout in {'uni-column', 'multi-column'}:
        extra_col_dict_list = []
    tt0 = time.time()
    tt1 = tt0
    newsflash('Processing input records ...')
    for in_tuple in panda_input.itertuples():
        """ Need to convert back to regular tuple, from pandafied named tuple with extra leading index number """
        in_supple = tuple(in_tuple[1:])
        # source_string = in_tuple[colname]
        source_string = in_supple[colno]
        source_terms = list(map(lambda x: x.strip(), source_string.split(",")))
        target_groups = {}
        if table_layout == 'in-situ' and keep_original:
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

        if table_layout in {'in-situ', 'uni-row', 'uni-column'}:
            target_string = ', '.join(tg_series.values)
            # newsflash("Additions to %s from MeSH: %s" % (source_string, target_string))
            # newsflash()

        if table_layout in {'in-situ', 'uni-row', 'multi-row'}:
            out_dict = dict(zip(out_columns, in_supple))
            # newsflash(out_dict)

            if table_layout == 'in-situ' or keep_original:
                if len(tg_series) > 0 or keep_original:
                    if table_layout == 'in-situ':
                        out_dict[colname] = target_string
                    out_dict_list.append(out_dict)

            if table_layout == 'uni-row' and len(tg_series) > 0:
                out_dict_extra = dict(out_dict)
                out_dict_extra[colname] = target_string
                out_dict_list.append(out_dict_extra)

            elif table_layout == 'multi-row':
                for hit in tg_series.values:
                    out_dict_iter = dict(out_dict)
                    out_dict_iter[colname] = hit
                    out_dict_list.append(out_dict_iter)

        elif table_layout in {'uni-column', 'multi-column'} and (len(tg_series) > 0 or keep_original):
            """ Now need to handle uni- and multi-column outputs """
            out_dict = dict(zip(out_columns, in_supple))
            out_dict_list.append(out_dict)
            extra_col_dict = {}
            if table_layout == 'multi-column':
                extra_col_dict = dict(tg_series)
            elif table_layout == 'uni-column':
                # extra_col_dict = dict({'all_ontologies', ', '.join(tg_series.values)})
                extra_col_dict = dict({'EQUIVALENT_TRAIT_URIS': target_string})
            extra_col_dict_list.append(extra_col_dict)

        in_tuple_counter += 1
        if in_tuple_counter % 4000 == 0:
            tt2 = time.time()
            newsflash("Processed %d thousand input records: took %.2f s (increment of %.2f s) ..." %
                           (int(in_tuple_counter / 1000), float(tt2 - tt0), float(tt2 - tt1)))
            tt1 = tt2

    tt2 = time.time()
    newsflash("Processed a total of %d input records, in %.2f seconds" % (in_tuple_counter, float(tt2 - tt0)))
    newsflash("No. of records in output spreadsheet is %d" % len(out_dict_list))
    newsflash()
    panda_output = pd.DataFrame(out_dict_list, columns=out_columns)
    if table_layout in {'uni-column', 'multi-column'}:
        newsflash("Adding new columns ...")
        # tg_series out of scope here!
        # extra_columns_df = pd.DataFrame(extra_col_dict_list,
        #                                 columns=[colname] if table_layout == 'uni-column' else tg_series.keys())
        extra_columns_df = pd.DataFrame(extra_col_dict_list)
        panda_output = pd.concat([panda_output.loc[:, :colname if keep_original else prev_colname], extra_columns_df,
                                  panda_output.loc[:, next_colname:]], axis=1)

    return panda_output


def re_ontologise(input_file, output, layout, file_format, column_index, keep, target, uri_format, boundary, paxo, oxo_url,
                  number, verbose):

    target = sorted(target)
    # newsflash("Length of target ontology array is %d" % len(target))
    # for t in target:
    #     newsflash("Target is %s" % t)
    field_separator = ',' if file_format == 'csv' else '\t'
    ss_dict = parse_ss(input_file, field_separator, column_index)
    iri_map = ss_dict['unique_iris']
    """ Print out list of source iris in iri_map """
    # iri_counter = 0
    # for src_iri in iri_map:
        ### print("%d\t%s\t%s" % (iri_counter, src_iri, iri_map[src_iri]))  # Print values _and_ keys
        # newsflash("%d\t%s" % (iri_counter, src_iri))
        # iri_counter += 1
    panda_original = ss_dict['pandafued']
    newsflash("Calling map_iris with url = '%s' ..." % oxo_url)
    map_iris(iri_map, target, boundary, paxo, oxo_url, number, verbose)
    newsflash("Calling augment ...")
    ontologically_enriched = augment(panda_original, iri_map, layout, column_index, keep, uri_format)
    """ Print out augmented_panda here ... """
    # newsflash("No. of dictionary elements: %d" % len(ss_dict))
    # newsflash("No. of rows in spreadsheet: %d" % len(panda_original))
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
    newsflash("Outputting ontologically enriched spreadsheet ...")
    # print(ontologically_enriched.head(30).to_csv(index=False, sep='\t'))
    print(ontologically_enriched.to_csv(index=False, sep=field_separator))


def main():

    """ First of all, check configuration file """
    parser1 = argparse.ArgumentParser(description='Config filepath', add_help=False)
    parser1.add_argument('-g', '--config', help='filepath of config file (ini format)')
    namespace, extra = parser1.parse_known_args()
    config_file = namespace.config

    """ Second of all, get configuration info from configuration file (if specified) """
    ontoconfig = configparser.ConfigParser(os.environ)
    ontoconfig.optionxform = str
    target_list = None
    if config_file is not None:
        ontoconfig.read(config_file)
        targets_plus = ontoconfig.items('Targets')
        spurious_targets = ontoconfig.items('DEFAULT')
        real_targets = dict(set(targets_plus) - set(spurious_targets))
        # newsflash(pd.Series(real_targets))
        target_list = list(real_targets.values())

    """ Third of all, parse the rest of the switches, possibly using defaults from configuration file """
    parser2 = argparse.ArgumentParser(prog='Ontomapper', description='Type of URI to output, etc.', parents=[parser1])
    cfg_sect_lookup = config_or_bust(ontoconfig, 'Params')
    parser2.add_argument('-i', '--input-file', default=cfg_sect_lookup('input_file', 'string'),
                         help='location of input spreadsheet: accepts filepath or URL')
    parser2.add_argument('-o', '--output', default=cfg_sect_lookup('output', 'string'),
                         help='output spreadsheet filepath **NO CURRENT EFFECT**')
    parser2.add_argument('-f', '--file-format', choices=['csv', 'tsv'],
                         default=cfg_sect_lookup('file_format', 'string'),
                         help='file format (both input and output)')
    parser2.add_argument('-l', '--layout', choices=['in-situ', 'uni-column', 'multi-column', 'uni-row', 'multi-row'],
                         default=cfg_sect_lookup('layout', 'string'),
                         help="%s%s" % ('whether new ontology terms are required in multiple rows, ',
                                        'multiple columns, a single row, a single column, or the originating cell'))
    parser2.add_argument('-c', '--column-index', type=int, default=cfg_sect_lookup('column_index', 'int'),
                         help='zero-based index of column containing source ontology terms')
    kmeg = parser2.add_mutually_exclusive_group(required=False)
    kmeg.add_argument('-k', '--keep', dest='keep', action='store_true', help='retain source ontology terms')
    kmeg.add_argument('-d', '--no-keep', dest='keep', action='store_false', help='ditch source ontology terms')
    # parser2.set_defaults(keep=cfg_sect_lookup('keep'))
    parser2.add_argument('-t', '--target', nargs='+', default=target_list,
                         help='space-separated list of target ontology prefixes')
    parser2.add_argument('-u', '--uri-format', choices=['long', 'short', 'curie'],
                         default=cfg_sect_lookup('uri_format', 'string'),
                         help='format of target ontology term identifiers **NO CURRENT EFFECT**')
    parser2.add_argument('-b', '--boundary', type=int, default=cfg_sect_lookup('boundary', 'int'),
                         help="%s%s" % ('minimum percentage confidence threshold of target ontology term matches ',
                                        '**NO CURRENT EFFECT: ENFORCE 100%% CONFIDENCE (OxO distance=1)**'))
    parser2.add_argument('-r', '--oxo-url', default=cfg_sect_lookup('oxo_url', 'string'),
                         help='OxO (or Paxo) web service URL')
    # parser2.add_argument('-p', '--paxo', type=bool, nargs='?', const=True, default=False,
    #                      help='Use Paxo web service rather than OxO **NO CURRENT EFFECT**')
    pmeg = parser2.add_mutually_exclusive_group(required=False)
    pmeg.add_argument('-p', '--paxo', dest='paxo', action='store_true', help='use Paxo rather than OxO')
    pmeg.add_argument('-x', '--no-paxo', dest='paxo', action='store_false', help='do not use Paxo: use OxO')
    parser2.add_argument('-n', '--number', type=int, default=cfg_sect_lookup('query_term_number', 'int'),
                         help='number of query terms to chunk, per HTTP request on the OxO web service')
    vmeg = parser2.add_mutually_exclusive_group(required=False)
    vmeg.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help="%s%s" % ('send verbose progess reports to standard error: ',
                                     'not recommended for regular use'))
    vmeg.add_argument('-q', '--quiet', dest='verbose', action='store_false', help='suppress verbose output')
    parser2.set_defaults(keep=cfg_sect_lookup('keep', 'boolean'), paxo=cfg_sect_lookup('paxo', 'boolean'),
                         verbose=cfg_sect_lookup('verbose', 'boolean'))
    parser2.add_argument('--version', action='version', version='%(prog)s 1.0')

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
    for inactive_arg in ['output', 'paxo', 'uri_format', 'boundary']:
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

    newsflash(arg_dict, arg_dict['verbose'])

    """ '**' unpacks a dictionary """
    re_ontologise(**arg_dict)


if __name__ == "__main__":
    main()
