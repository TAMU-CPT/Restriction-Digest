from bs4 import BeautifulSoup
import sys
with open(sys.argv[1],'r') as file:
	soup = BeautifulSoup(file)
	list2 = soup.find_all('i')
	for j in range(30):
		print j,list2[j]
	list1  =  soup.find_all('code')
	for i in range(10):
		print i,list1[i]
