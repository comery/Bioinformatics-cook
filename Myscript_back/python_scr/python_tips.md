### pirnt infomation
##### [time](http://www.runoob.com/python/python-date-time.html)
```python
#!/usr/bin/python
# -*- coding: UTF-8 -*-
 
import time
 
# 格式化成2016-03-20 11:45:39形式
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
 
# 格式化成Sat Mar 28 22:24:24 2016形式
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()) )
  
# 将格式字符串转换为时间戳
a = "Sat Mar 28 22:24:24 2016"
print(time.mktime(time.strptime(a,"%a %b %d %H:%M:%S %Y")))

```
output:
```text
2016-04-07 10:25:09
Thu Apr 07 10:25:09 2016
1459175064.0
```

##### Calendar
```python
#!/usr/bin/python
# -*- coding: UTF-8 -*-
 
import calendar
 
cal = calendar.month(2016, 1)
print "以下输出2016年1月份的日历:"
print cal
```
output:
```text
以下输出2016年1月份的日历:
    January 2016
Mo Tu We Th Fr Sa Su
             1  2  3
 4  5  6  7  8  9 10
11 12 13 14 15 16 17
18 19 20 21 22 23 24
25 26 27 28 29 30 31
```

### control process

#### os._exit() & sys.exit() & exit() & exit(1)

```python
#!/usr/local/bin/env python
import os, sys

try:
    sys.exit(0)
except:
    print('die')
finally:
    print('cleanup')

try:
    os._exit(0)
except:
    print('die')
print('os.exit')#不打印直接退出了
```

output:
```text
die
cleanup
```
>区别

>综上，sys.exit()的退出比较优雅，调用后会引发SystemExit异常，可以捕获此异常做清理工作。os._exit()直接将python解释器退出，余下的语句不会执行。

>一般情况下使用sys.exit()即可，一般在fork出来的子进程中使用os._exit()
一般来说os._exit() 用于在线程中退出 
sys.exit() 用于在主线程中退出。
exit() 跟 C 语言等其他语言的 exit() 应该是一样的。 
os._exit() 调用 C 语言的 _exit() 函数。
builtin.exit 是一个 Quitter 对象，这个对象的 call 方法会抛出一个 SystemExit 异常。

>exit(0)和exit(1)
exit(0)：无错误退出 
exit(1)：有错误退出 
退出代码是告诉解释器的（或操作系统）

#### subprocess.check_output() 获取外部命令的结果
```python
less_cmd = "less -S " + filename + "|wc -l"
priwc = subprocess.check_output(less_cmd,shell=True)
primer_lines = priwc.decode('utf-8')
```

#### 反向互补序列
```python
def comp_rev(sequence):
    # make a sequence complement #
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()[::-1]
    
def reverseComplement(s):
    complement = { 'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A' }
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def reverseComplement(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence
```

​        

#### 二维字典
建立一个二维字典

```python
dict_2d = {'a': {'a': 1, 'b': 3}, 'b': {'a': 6}}

```
访问这个字典的元素
```python
dict_2d['a']['b']
```
但是”2-D” dictionary 新添一个”key-value”对时，不能简单的用的形式。
```python
dict_2d['a']['c'] = 8
```
因为二维字典的两层key和value之间会混淆，需要判断第一个key是否已经存在了。添加二维的字典可以通过一个函数来简单实现：
```python
def addtwodimdict(thedict, key_a, key_b, val): 
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})
```



#### 指定运行python版本

>如果你用 python xxoo.py 来运行，那么写不写都没关系，如果要用 ./xxoo.py 那么就必须加这行，这行被称为 shebang, 用来为脚本语言指定解释器.通常认为用 

```python
#!/usr/bin/env python 
```
要比 
```python
#!/usr/bin/python 
```
更好，因为 python 解释器有时并不安装在默认路径，例如在 virtualenv 中。


#### read fasta
```python
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, "''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, "''.join(seq))

with open('f.fasta') as fp:
    for name, seq in read_fasta(fp):
        print(name, seq)
```

```python
import pyfastx
if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} *.fasta")

fa = pyfastx.Fastx('*.fasta')
for name,seq,comment in fa:l

    if name in pair:
        print(f">{pair[name]}\n{seq}")
```

#### read fastq
```python
def parse_pe_fastq(fq1, fq2):
    while True:
        name1 = fq1.readline().split(" ")[0]
        name2 = fq2.readline().split(" ")[0]
        if not name1 or not name2:
            break
        read1, nothing1, qual1 = fq1.readline()[:-1], fq1.readline(), fq1.readline()[:-1]
        read2, nothing2, qual2 = fq2.readline()[:-1], fq2.readline(), fq2.readline()[:-1]
        assert name1 == name2, 'fastq1, fastq2 is not paired-end'
        yield name1, read1, read2, qual1, qual2
```


#### open glob files
```python
import glob,os
 
if __name__=='__main__':
    prefix=input('Input the prefix of images:')
    files=glob.glob(prefix+'_*')
    num=len(files)
```


## convert phylip to fasta
```python
def phylip2fasta(inpath, outfile):
    try:
        infile = open(inpath,"r")
    except OSError:
        exit("There was a problem opening the specified directory. Are you sure it exists? Quitting.")

    out = open(outfile,"w")

    firstline = infile.readline()
    while True:
        name = infile.readline().strip()
        seq = infile.readline().strip()
        if name:
            out.write(">" + name + "\n")
            out.write(seq + "\n")
        else:
            break

    infile.close()
    out.close()
```



### read fastq

```python
import mappy as mp
with open('out.fq', 'w') as out:
	for name, seq, qual, comment in mp.fastx_read('./test.fq', read_comment=True):
		out.write(f"{read[0]}")
```



