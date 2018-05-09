#!/usr/bin/env python3

import sys


def newsflash(msg=None):
    """
    Sends a message to the standard error console; gets around posting info
    to console while redirecting standard output to a file (or whatever).
    """
    if msg is None:
        msg = ""
    sys.stderr.write("%s\n" % (msg))


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
