fa_in = open("in.fa", "r")　　#读入FASTA文件
fa_Info = []　　#定义一个列表，用于存储序列的信息
fa_Seq = []　　#定义一个列表，用于存储序列
fa_Num = -1
for line in fa_in.readlines():　　#每次从FASTA文件中读取一行
    line = line.rstrip()　　#去掉行末的换行符
    if line[0] == ">":　　#判断如果是信息行，就存储在fa_Info，否则存储在fa_Seq
        fa_Info.append(line)
        fa_Num = fa_Num + 1
        fa_Seq.append("")
    else:
        fa_Seq[fa_Num] = fa_Seq[fa_Num] + line
 
fa_out = open("out.fa", "w")　　#打开一个文件，用于写入
for i in range(fa_Num + 1):
    fa_out.write(fa_Info[i] + "\n")
    while len(fa_Seq[i]) > 60:　　#每行写60个碱基
        fa_out.write(fa_Seq[i][:60] + "\n")
        fa_Seq[i] = fa_Seq[i][60:]
    else:
        fa_out.write(fa_Seq[i] + "\n")
