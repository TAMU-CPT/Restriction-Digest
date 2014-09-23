from bs4 import BeautifulSoup
import sys
with open(sys.argv[1],'r') as file:
	soup = BeautifulSoup(file)
	list = soup.find_all('tr')
	for element in list:
		Temp_Soup = BeautifulSoup(str(element))
		Temp_List = Temp_Soup.find_all('td')
		for e in Temp_List:
			Sub_Soup = BeautifulSoup(str(e))
			print Sub_Soup.get_text()
	
	
