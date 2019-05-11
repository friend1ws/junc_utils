#! /usr/bin/env python3

# from parser import create_parser
# from . import parser.create_parser
# import parser as ppp
from . import parser

def main():

    cparser = parser.create_parser()
    args = cparser.parse_args()
    args.func(args)

