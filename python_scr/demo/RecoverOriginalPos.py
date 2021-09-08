#!/usr/bin/env python3
import sys
import os

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



def recovered_pos(agp_file):
    original_pos = {}
    first_region = {}
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
                # map new(chr) location into old(scaffold) location
                if strand == '+':
                    for i in range(super_scaf_start-1, super_scaf_end):
                        real_pos = i + 1
                        real_block_pos = real_pos - super_scaf_start + block_start
                        addtwodimdict(original_pos, super_scaf, real_pos, (block_id, real_block_pos))
                else:
                    for i in range(super_scaf_start-1, super_scaf_end):
                        real_pos = i + 1
                        y1 = real_pos - super_scaf_start + block_start # position by block
                        y0 = block_end - block_start + 2 - y1  # reverse complement
                        addtwodimdict(original_pos, super_scaf, real_pos, (block_id, y0))
            else:
                sys.exit("unrecognized type")

    return original_pos


def main():
    original_pos = recovered_pos(sys.argv[1])
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            ref_id = tmp[0]
            ref_start = int(tmp[1])
            ref_end = int(tmp[2])
            que_id = tmp[3]
            que_start = int(tmp[4])
            que_end = int(tmp[5])
            output = []
            # if this chr is makeuped, then retrieve original position
            if ref_id in original_pos.keys():
                original_ref_s = original_pos[ref_id][ref_start][0]
                original_ref_e = original_pos[ref_id][ref_end][0]
                if original_ref_s != original_ref_e:
                    print("this sv span over two blocks, remove")
                    continue
                else:
                    original_ref_start = original_pos[ref_id][ref_start][1]
                    original_ref_end = original_pos[ref_id][ref_end][1]
                    output = [original_ref_s, original_ref_start, original_ref_end]
            # if this chr is not makeuped, do nothing
            else:
                output = tmp[:3]

            if que_id in original_pos.keys():
                original_que_s = original_pos[que_id][que_start][0]
                original_que_e = original_pos[que_id][que_end][0]
                if original_que_s != original_que_e:
                    print("this sv span over two blocks, remove")
                    continue
                else:
                    original_que_start = original_pos[que_id][que_start][1]
                    original_que_end = original_pos[que_id][que_end][1]
                    output.extend([original_que_s, original_que_start, original_que_end])
            else:
                output.extend(tmp[4:7])

            output.extend(tmp[6:])
            if output[1] > output[2]:
                output[1], output[2] = output[2], output[1]
            if output[4] > output[5]:
                output[4], output[5] = output[5], output[4]
            output = [str(n) for n in output]
            print("\t".join(output))


if __name__ == '__main__':
    main()
