#!/usr/env/python python3
import sys
import gzip
import mappy as mp

if len(sys.argv) != 4:
    usage = """
    Downsample fastq file by size(bp/kb/mb/gb) or line or reads,
    you can specify the size you want to generate, like 1g, 1m, 1k,
    or 10000l, which mean 1 gigabase, 1 million base, 1 kilobase, or
    10000 line. if your data is pair-end, you'd better deal R1,R2 together
    """
    print(usage)
    print("Usage: python3 {} *_R1.fastq.gz 10000l out.fq.gz".format(sys.argv[0]))
    print("Usage: python3 {} *_R1.fastq.gz,*_R2.fastq.gz 5g out_R1.fq.gz,out_R2.fq.gz".format(sys.argv[0]))
    exit()


def get_size(s):
    if s.endswith("g"):
        size = float(s.replace("g", "")) * 10 ** 9
    elif s.endswith("m"):
        size = float(s.replace("m", "")) * 10 ** 6
    elif s.endswith("k"):
        size = float(s.replace("k", "")) * 10 ** 3
    elif "." in s or '-' in s:
        sys.exit(f"can not recognize size type of {s}")
    else:
        size = int(s)
    return size

def smart_open(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')


def main():
    input_fq = sys.argv[1]
    if "," in input_fq:
        inputs = input_fq.split(",")
        if inputs[0] == inputs[1]:
            sys.exit(f"how careless you are! two input fastqs are same!")

        if len(inputs) == 2 and inputs[1] != "":
            mode = 'pe'
        else:
            mode = 'se'
    else:
        mode = 'se'

    if mode == 'se': # single-end mode
        out = gzip.open(sys.argv[3], 'wb')
        if sys.argv[2].endswith("l"):
            lines_required = int(sys.argv[2].replace("l", ""))
            #print(lines_required)
            # sampling by lines
            current_line = 0
            current_base = 0
            for read in mp.fastx_read(input_fq, read_comment=False):
                current_line += 1
                current_base += len(read[1])
                if current_line <= lines_required:
                    tmp = f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n"
                    out.write(tmp.encode())
                else:
                    real_line = current_line - 1    # roll back
                    real_base = current_base - len(read[1])
                    print(f"output generated: {real_line}\t{real_base}")
                    break
        else:
            size = get_size(sys.argv[2])
            #print(size)
            current_size = 0
            for read in mp.fastx_read(input_fq, read_comment=False):
                current_size += len(read[1])
                if current_size <= size:
                    tmp = f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n"
                    out.write(tmp.encode())
                else:
                    real_size = current_size - len(read[1]) # roll back
                    print(f"output generated: {real_size}")
                    break

        out.close()
    else: # pair-end mode
        outs= sys.argv[3].split(",")
        if len(outs) != 2 or outs[1] == "":
           sys.exit(f"Given you give two input fastq files, you must give two output files")

        out1 = gzip.open(outs[0], 'wb')
        out2 = gzip.open(outs[1], 'wb')

        if sys.argv[2].endswith("l"):
            lines_required = int(sys.argv[2].replace("l", ""))
            #print(lines_required)
            # sampling by lines
            current_line = 0
            current_base = 0
            fq1 = mp.fastx_read(inputs[0], read_comment=False)
            fq2 = mp.fastx_read(inputs[1], read_comment=False)
            for read1, read2 in zip(fq1, fq2):
                current_line += 1
                current_base += len(read1[1]) + len(read2[1])
                if current_line <= lines_required:
                    tmp1 = f"@{read1[0]}\n{read1[1]}\n+\n{read1[2]}\n"
                    tmp2 = f"@{read2[0]}\n{read2[1]}\n+\n{read2[2]}\n"
                    out1.write(tmp1.encode())
                    out2.write(tmp2.encode())
                else:
                    # roll back
                    real_line = current_line - 1
                    real_base = current_base - len(read1[1]) - len(read2[1])
                    print(f"output generated: {real_line}\t{real_base}")
                    break
        else:
            size = get_size(sys.argv[2])
            #print(size)
            current_size = 0
            fq1 = mp.fastx_read(inputs[0], read_comment=False)
            fq2 = mp.fastx_read(inputs[1], read_comment=False)
            for read1, read2 in zip(fq1, fq2):
                current_size += len(read1[1]) + len(read2[1])
                if current_size <= size:
                    tmp1 = f"@{read1[0]}\n{read1[1]}\n+\n{read1[2]}\n"
                    tmp2 = f"@{read2[0]}\n{read2[1]}\n+\n{read2[2]}\n"
                    out1.write(tmp1.encode())
                    out2.write(tmp2.encode())
                else:
                    real_size = current_size - len(read1[1]) - len(read2[1])
                    print(f"output generated: {real_size}")
                    break

        out1.close()
        out2.close()

if __name__ == '__main__':
    main()
