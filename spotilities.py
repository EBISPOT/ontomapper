#!/usr/bin/env python3

import sys
import configparser


def newsflash(msg=None, verbose=True):
    """
    Sends a message to the standard error console; gets around posting info
    to console while redirecting standard output to a file (or whatever).
    """
    if verbose:
        if msg is None:
            msg = ""
        sys.stderr.write("%s\n" % (msg))


def config_or_bust(cfg_object, cfg_section):

    def config_section_lookup(cfg_key, cfg_type):
        try:
            if cfg_type == 'string':
                cfg_value = cfg_object.get(cfg_section, cfg_key)
            elif cfg_type == 'int':
                cfg_value = cfg_object.getint(cfg_section, cfg_key)
            elif cfg_type == 'float':
                cfg_value = cfg_object.getfloat(cfg_section, cfg_key)
            elif cfg_type == 'boolean':
                cfg_value = cfg_object.getboolean(cfg_section, cfg_key)
            else:
                raise ValueError('cfg_type should be one of: string, int, float, boolean')
        except configparser.NoSectionError:
            cfg_value = None
        return cfg_value
    return config_section_lookup


def listify_uris(uri_string):
    """
    Take a comma-separated string of URIs and return a list of the same
    :param uri_string:
    :return uri_list:
    """


def stringify_uris(uri_list):
    """
    Take a list of URIs and return a comma-separated list of the same
    :param uri_list:
    :return uri_string:
    """


def extend_uri_list():
    """
    Add URI
    :param short_list long_list:
    :return count_of_records_added_to_long_list:
    """
