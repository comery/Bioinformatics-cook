#!/usr/bin/python3
f=open("test.fa",'r')
#print (f.tell())
a = f.readline(4)
print(a)
#print (f.tell())
for line in f:
    line = line.strip()
#    print(line," ",len(line),end='\n')

#try:
#    fock = open ("1",'r')
#except :
#    print ("The file does not exist, exiting gracefully")

