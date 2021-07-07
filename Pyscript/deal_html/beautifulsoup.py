from bs4 import BeautifulSoup
import lxml
import sys

def delete(str):
	    arr = str.split("\n")
	    new = ""
	    for i in arr:
	        if not len(i) == 0:
	            new += i +"\n"
	    return new

if len(sys.argv) < 2:
	print("Usage:\n\t python " + sys.argv[0] + " imput.html")
	exit()

with open(sys.argv[1],'r') as fh:
	html = fh.read()
	soup = BeautifulSoup(html,"lxml")
	'''print(soup.prettify())'''
	print(soup.p.string)
	for pp in  soup.find_all('p'):
		print(pp.string)
	'''
	print("before removing nona line:\n")
	print(soup.get_text())
	print("after removing nona line:\n")
	print(delete(soup.get_text()))
        '''
