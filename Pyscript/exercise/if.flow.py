#!/usr/bin/python
x= int(input("Please input an integer: " ))
if x < 0:
    x = 0
    print('Negative change to Zero')
elif x == 0:
    print('Zero')
elif x == 1:
    print('Single')
else:
    print('More')


