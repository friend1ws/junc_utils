#! /usr/bin/env python

from . import parser

def main():

    cparser = create_parser()
    args = cparser.parse_args()

    if vars(args) == {}:
        parser.print_usage()
    else:
        args.func(args)

