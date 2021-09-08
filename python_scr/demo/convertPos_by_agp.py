#!/usr/bin/env python3
import sys
import os
#from icecream import ic

desc = """
the var file looks like:
ref_id  ref_start   ref_end qry_id  qry_start   qry_end
"""
if len(sys.argv) < 3:
    print(desc)
    sys.exit("python3 {}  <*.agp> <*.var>".format(sys.argv[0]))


def addtwodimdict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})


def function_by_region(agp_file):
    rules = {}
    coords = {}
    with open(agp_file, 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            super_scaf = tmp[0]
            super_scaf_start = int(tmp[1])
            super_scaf_end = int(tmp[2])
            order = int(tmp[3])
            block_type = tmp[4]
            block_id = tmp[5]
            block_start = int(tmp[6])
            block_end = int(tmp[7])
            strand  = tmp[8]
            if block_type == 'N':
                continue
            elif block_type == 'W':
                # key_a is current scaffold id
                # key_b is the region position (start, end)
                # val is coording block information
                addtwodimdict(coords, super_scaf, (super_scaf_start, super_scaf_end), (block_id, block_start, block_end, strand))
                assert abs(super_scaf_end-super_scaf_start) == abs(block_end-block_start), "bad length"
                # map new(chr) location into old(scaffold) location, generate a rule
                def rule(rs, re, bid, bs, be, strand, x):
                    """
                    x is one position that you want to convert
                    """
                    phase = x - rs
                    if strand == '+':
                        target = bs + phase
                    else:
                        target = be - phase
                    return bid, target
                # key_a is current scaffold id
                # key_b is the region position (start, end)
                # val is a function
                addtwodimdict(rules, super_scaf, (super_scaf_start, super_scaf_end), rule)

            else:
                sys.exit("unrecognized type")

    return rules, coords


def which_fuction(rules, coords, super_scaf, search):
    """
    which region does a position locate in ?
    """
    sign = 0
    end = 0
    for i in rules[super_scaf].keys():
        # block start and end
        (s, e) = i
        (bid, bs, be, strand) = coords[super_scaf][i]
        end = e
        if search >= s and search <= e: # in this block, so found it
            sign = 1
            break
    if sign:
        return s, e, bid, bs, be, strand, rules[super_scaf][i]
    elif search > end:
        print(f"{super_scaf} {search} beyond the entire block")
        return None, None, None, None, None, None, None
    else:
        # N regions
        print(f"can not locate this pos {super_scaf} {search}, it might be in gap")
        return None, None, None, None, None, None, None


def main():
    rules, coords = function_by_region(sys.argv[1])
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            ref_id = tmp[0]
            if ref_id not in coords.keys():
                # this scaffold was not makeuped, no change
                print(f"{ref_id}\t{tmp[1]}")
            else:
                search = int(tmp[1])
                rs, re, bid, bs, be, strand, function = which_fuction(rules, coords, ref_id, search)
                if function:
                    original_scaf, original_pos = function(rs, re, bid, bs, be, strand, search)
                    print(f"{original_scaf}\t{original_pos}")


if __name__ == '__main__':
    main()
