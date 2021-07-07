#!/usr/bin/env python3
import sys
if len(sys.argv) < 5:
    print("Usage: python3 {} {}".format(sys.argv[0], "<source> <field> <specific item> <field>"))
    exit()

def get_list(File, field):
    records = {}
    with open(File, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            records[tmp[field-1]] = 1

    return records


def main():
    records = get_list(sys.argv[3], int(sys.argv[4]))
    field_target = int(sys.argv[2])
    found = []
    err = open(sys.argv[3] + ".missing.out", 'w')
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            if tmp[field_target-1] in records.keys():
                print(i.strip())
                found.append(tmp[field_target-1])
    for k in records.keys():
        if k not in found:
            err.write(k + "\n")
    err.close()


if __name__ == '__main__':
    if sys.argv[2] == '0' or sys.argv[4] == '0':
        print("field is 1-based")
        exit()
    main()
