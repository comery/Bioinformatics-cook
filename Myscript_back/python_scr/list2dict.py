from time import time
t = time()

list = ["hello", "my", "python", "world"]

list = dict.fromkeys(list, True)

filter = []
for i in range(100000):
    for find in ['a', 'is', 'b', 'my']:
        if find in list:
           filter.append(find)

print(time() - t)
