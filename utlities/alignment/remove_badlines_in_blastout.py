#!/usr/bin/evn python3

import sys
#import argsparse
des = "This script is used in a specific condition. when you run blast,\n" +\
        "and it stopped by accident, now you want to run remaining sequence,\n" +\
        "now you need to split them out, and run again. finally you need merge\n" +\
        " them together, in that time, you need this script to deal bad lines.\n"
if len(sys.argv) < 3:
    print(des)
    print("Usage: python3 {} {} {}".format(sys.argv[0],
                                              "<blastout|file>",
                                              "<output>"))
    exit()

def parser_blast(file):
    pre = ""
    target_match = []
    with open(file, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if not line:
                break
            tmp = line.strip().split("\t")
            if pre == "":
                target_match = [line,]
                pre = tmp[0]
            elif tmp[0] == pre:
                target_match.append(line)
            else:
                yield target_match
                target_match = [line,]
                pre = tmp[0]
        yield target_match


def main():
    with open(sys.argv[2], 'w') as fw:
        for arr in parser_blast(sys.argv[1]):
            record = []
            for line in arr:
                tmp = line.split("\t")
                if len(tmp) < 12:
                    continue
                if line not in record:
                    print(line, file=fw)
                    record.append(line)

if __name__ == "__main__":
    main()


