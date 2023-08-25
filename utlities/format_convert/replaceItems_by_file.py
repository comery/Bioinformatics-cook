#!/usr/bin/env python3
import sys
if len(sys.argv) < 4:
    help="""
    convert specific fields in your target file by a list which contains cooridinates
    """
    print(help)
    print("Usage: python3 {} {} {} {}".format(sys.argv[0], "target_file", "columnN", "convert_list"))
    exit()

def get_list(File):
    records = {}
    with open(File, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            records[tmp[0]] = tmp[1]

    return records


def main():
    records = get_list(sys.argv[3])
    field_target = int(sys.argv[2])
    found = []
    err = open(sys.argv[3] + ".missing.out", 'w')
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            if field_target > len(tmp):
                sys.exit(f"your target file only has {len(tmp)} fields! Can not find field [{field_target}]")
            if tmp[field_target-1] in records.keys():
                tmp[field_target-1] = records[tmp[field_target-1]]
                print("\t".join(tmp))
            else:
                print("# it will not be conveted:", file=err)
                print(i.strip(), file=err)


if __name__ == '__main__':
    if sys.argv[2] == '0':
        print("field is 1-based")
        exit()
    main()
