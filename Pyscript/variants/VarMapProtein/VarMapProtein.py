#!/usr/bin/env python3

import sys
import argparse


class Feature:

    def __init__(self, )













def main(args):

    with open(args.link, 'r') as fh:


if __name__ == '__main__':
    des = """

    xxx

    yangchentao at genomics.cn, BGI.
    """
    parser = argparse.ArgumentParser(description=des,
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-g", "--genome", type=str, required=True, metavar='<STR>',
                        help="genome seuqnce")
    parser.add_argument("-f", "--feature", type=str, required=True, metavar='<STR>',
                        help="gene feature, gff/gtf/bed format")
    parser.add_argument("-v", "--var", type=str, required=True, metavar='<STR>',
                        help="variants files")
    parser.add_argument("-l", "--long", action="store_true",
                        help="long information")
    args =  parser.parse_args()
    main(args)
