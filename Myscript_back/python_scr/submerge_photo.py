import matplotlib.pyplot as plt

plt.figure()

from PIL import Image

arr = ['p1.png', 'p2.png', 'p3.png', 'p4.png']
toImage = Image.new('RGBA',(400,400))
for i in range(4):
    fromImge = Image.open(arr[i])
    # loc = ((i % 2) * 200, (int(i/2) * 200))
    loc = ((int(i/2) * 200), (i % 2) * 200)
    print(loc)
    toImage.paste(fromImge, loc)

toImage.save('merged.png')

