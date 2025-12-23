#-*- coding:utf-8 -*-
import sys
import os

if __name__ == '__main__':

    #sys.argv[0]这个参数是脚本名称.
    if len(sys.argv) < 3:
        sys.exit("usage: python xxx.py sourceFile  outFile  .....")

    sourceFilePath = sys.argv[1];   #参数1
    outFilePath    = sys.argv[2];   #参数2
    if os.path.exists(sourceFilePath):

        with open(sourceFilePath , "r") as f:
            sourceFileStr = f.read();
        #四个空格
        outFileStr = sourceFileStr.replace("\t","    ")

        with open(outFilePath,"a+") as w:
            w.write(outFileStr);
        print("replace success! , output path : " + outFilePath)
    else:
        print("usage: python xxx.py sourceFile  outFile  .....")
