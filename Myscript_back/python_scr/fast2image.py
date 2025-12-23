import sys
import cv2
from PIL import Image
import numpy as np
from scipy import misc

image = Image.open("test.png")  # 首先在该py文件所在目录下随便放一张图片，使用PIL.Image库的open方法打开
image_array = np.array(image)  # 使用numpy将该图片的二进制数据转换成多维数组形式
print(image_array)
misc.imsave('out.jpg', image_array)  # 使用misc.imsave方法将数组保存为图片

# make a dictory of atcg string with number (0-255)
base_code = {}
bases = ['A', 'T', 'C', 'G']
code = 0
for a in bases:
    for b in bases:
        for c in bases:
            for d in bases:
                string = a + b + c + d
                base_code[string] = code
                code += 1



def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


with open(sys.argv[1] ,'r') as fp:
    for name, seq in read_fasta(fp):
        seq = seq.replace("\n", "")
        tmp = 0
        length= len(seq)
        new = ""
        pix = []
        R = []
        G = []
        B = []
        while tmp < length:
            new = seq[tmp:tmp+12]
            if len(new) == 12:
                tmp += 12
                str_b = new[0:4]
                str_g = new[4:8]
                str_r = new[8:12]
                dot = [base_code[str_b], base_code[str_r], base_code[str_r]]
                B.append(base_code[str_b])
                R.append(base_code[str_r])
                G.append(base_code[str_g])
                pix.append(dot)
            else:
                break
        print(len(pix))
        image_array = np.array(pix)
        image_array = image_array.reshape(32,12)
        misc.imsave('out.jpg', image_array)
        cv2.imwrite('out1.png', image_array)
        #R、G、B分量的提取
        #(B,G,R) = cv2.split(image)#提取R、G、B分量
        #R、G、B的合并
        merged = cv2.merge([np.array(B), np.array(G), np.array(R)])#合并R、G、B分量
        print(merged)
        cv2.imwrite("out2.png",merged)

