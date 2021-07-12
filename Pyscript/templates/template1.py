#!/usr/bin/env python3

import sys
import argparse

# Import palettable
try:
    from palettable.tableau import GreenOrange_12
except ImportError:
    print("We need palettable, sorry...\ninstall: pip install palettable")
    sys.exit(1)

try:
    import palettable.cartocolors.qualitative
except ImportError:
    print("We need palettable, sorry...\ninstall: pip install palettable")
    sys.exit(1)


def main(args):

    with open(args.link, 'r') as fh:

if __name__ == '__main__':
    des = """

    xxx

    yangchentao at genomics.cn, BGI.
    """
    parser = argparse.ArgumentParser(description=des,
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--link", type=str, required=True, metavar='<STR>',
                        help="BASTA annotation file(s) separated by comma")
    parser.add_argument("--scaf_len", type=str, required=True, metavar='<STR>',
                        help="all scaffold length, including querry and ref")
    parser.add_argument("--minL", type=int, required=False, metavar='<INT>',
                        help="min scaffold lenth in link record, default=1000000",
                        default=1000)
    parser.add_argument("--minB", type=int, required=False, metavar='<iNT>',
                        help="min block length of synteny, default=100000",
                       default=500)
    parser.add_argument("--rate", type=float, required=False,
                        metavar='<FLOAT>',
                        help="remove these scafoold which synteny block is shorter than scaf_len * rate",
                       default=0)
    args =  parser.parse_args()
    main(args)
