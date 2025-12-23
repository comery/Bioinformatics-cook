import re
 
#读取基因组序列，将碱基信息存储到字典中
genome_dict = {}
genome = open('Bacillus_subtilis.scaffolds.fasta', 'r')
for line in genome:
    line = line.strip()
    if line[0] == '>':
        seq = line.split('>')[1]
        genome_dict[seq] = ''
    else:
        genome_dict[seq] += line
 
genome.close()
 
#读取滑窗位置信息，并统计每段滑窗的 GC 含量并追加在滑窗测序深度统计文件后方
output = open('depth_gc.txt', 'w')
print('seq\tstart\tend\tDepth\tGC', file = output)
 
depth = open('depth.txt', 'r')
for line in depth:
    line = line.strip().split('\t')
    GC = len(re.findall('[GCgc]', genome_dict[line[0]][int(line[1]):int(line[2])]))  (int(line[2])-int(line[1]))
    print(f'{line[0]}\t{int(line[1])+1}\t{line[2]}\t{line[3]}\t{GC}', file = output)
 
depth.close()
output.close()