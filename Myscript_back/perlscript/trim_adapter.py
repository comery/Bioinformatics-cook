import sys
import re
from Bio import SeqIO

def rmPE(read1,read2,adaptor1,adaptor2,min_length):
    res_1 = rmSE(read1,adaptor1,min_length)
    res_2 = rmSE(read2,adaptor2,min_length)
    if res_1 and res_2:
        return res_1,res_2
    else:
        return False

def rmSE(read,adaptor,min_length):
    seq = read.seq
    seed_len = 6
    a_len = len(adaptor)
    seq_len = len(seq)
    for i in range(a_len - seed_len):
        seed = adaptor[i:i+seed_len]
        pos = 0
        while(pos < seq_len):
            find_pos = seq.find(seed,pos)
            if find_pos > 0:
                mistaken_count = 0
                _b = find_pos
                _e = find_pos + seed_len
                while(_b >= 0 and i >= find_pos - _b):
                    if adaptor[i - find_pos + _b] != seq[_b]:
                        mistaken_count += 1
                    if mistaken_count > 3:
                        break
                    _b -= 1
                else :
                    while(_e < seq_len and i - find_pos + _e < a_len):
                        if adaptor[ i - find_pos + _e ] != seq[_e]:
                            mistaken_count += 1
                        if mistaken_count > 3:
                            break
                        _e += 1
                    else:
                        if find_pos - i > min_length:
                            return  read[:find_pos-i]
                        else :
                            return False
                pos = find_pos + 1
            else:
                break
    return read

def rmAdaptor(argv):
    argv.pop(0)
    type = argv.pop(0)
    if type=='PE':
        read1_file,read2_file,adaptor1,adaptor2,out_prefix,min_length = argv
        read2_records = SeqIO.parse(open(read2_file),'fastq')
        read1_out = open( '%s.1.fq'%out_prefix,'w' )
        read2_out = open( '%s.2.fq'%out_prefix,'w' )
        for read1 in SeqIO.parse(open(read1_file),'fastq'):
            read2 = read2_records.next()
            rmPE_res = rmPE(read1,read2,adaptor1,adaptor2,min_length)
            if rmPE_res:
                read1_out.write(rmPE_res[0].format('fastq'))
                read2_out.write(rmPE_res[1].format('fastq'))
    elif type=='SE':
        reads_file,adaptor,out_prefix,min_length = argv
        reads_out = open( '%s.single.fq'%out_prefix,'w' )
        for reads in SeqIO.parse(open(reads_file),'fastq'):
            rmSE_res = False
            if re.search('[\s\/](\d)',reads.id).group(1) == '1':
                rmSE_res = rmSE(reads,adaptor1,min_length)
            elif re.search('[\s\/](\d)',reads.id).group(1) == '2':
                rmSE_res = rmSE(reads,adaptor2,min_length)
            if rmSE_res:
                reads_out.write(rmSE_res.format('fastq'))

if __name__ == '__main__':
    rmAdaptor(sys.argv)