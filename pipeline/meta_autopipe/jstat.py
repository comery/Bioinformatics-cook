#!/usr/bin/env python3
import json
import sys

if len(sys.argv) < 2:
    print("Usage: python3 {} {}".format(sys.argv[0], "xxx.json"))
    exit()


def load_from_db(handle):
    data = json.load(handle)
    return data

def main():
    with open(sys.argv[1], 'r') as fh:
        data = load_from_db(fh)
        work_status = data['work_status']
        print("*---Power by yangchentao---*")
        print("Status:\t2 = done, 1 = submitted, 0 = waitting\n")
        print("#ID\tfilter\tmerge\trmcondam\tassembly\tmapping")
        for i in work_status.keys():
            stat = [ str(x) for x in work_status[i] ]
            print(i + "\t" + "\t".join(stat))

if __name__ == "__main__":
    main()
