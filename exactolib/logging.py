#!/usr/bin/python3

"""
The purpose of this python3 script is to refine variants.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


import logging


def get_logger(name):
    logging.basicConfig(
        format = '%(asctime)s %(levelname)-8s %(message)s',
        level = logging.INFO,
        datefmt = '%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(name)

