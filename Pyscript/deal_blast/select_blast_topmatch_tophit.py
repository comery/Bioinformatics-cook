#!/usr/bin/evn python3

import sys
#import argsparse

if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {}".format(sys.argv[0],
                                              "<blastout|file>",
                                              "<topmatch|int>",
                                              "<tohit|int>"))
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
    match_lim = int(sys.argv[2])
    hit_lim = int(sys.argv[3])

    for arr in parser_blast(sys.argv[1]):
       # print(len(arr))
        match = 0
        hit = {}
        pre = " "
        output = []
        for line in arr:
            tmp = line.split("\t")
            querry = tmp[1]
            if pre == " ":
                match = 1
                pre = querry
                hit[pre] = 1
                output.append(line)
            elif querry == pre and hit[querry] >= hit_lim:
                continue
            elif querry == pre:
                hit[querry] += 1
                output.append(line)
            elif match >= match_lim:
                break
            else:
                match += 1
                pre = querry
                hit[pre] = 1
                output.append(line)

        print("\n".join(output))

if __name__ == "__main__":
    main()


