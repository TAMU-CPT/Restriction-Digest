from bs4 import BeautifulSoup
import sys
with open(sys.argv[1],'r') as file:

#    soup = BeautifulSoup(file)
#    list = soup.find_all('td')
#    new_list = []
#    new_new_list = []
#    for element in list:
       # Temp_Soup = BeautifulSoup(str(element))
       # Temp_List = Temp_Soup.find_all('td')
       # for e in Temp_List:
           # new_list += [BeautifulSoup(str(e)).get_text()]
    #Desired Table Headers, which should be the same on all the pages because they have the same format
#        new_list += [BeautifulSoup(str(element)).get_text()]
#    print new_list[17:23]
#    inc = 0
#    while 29+ inc <=len(new_list):
#        print new_list[23+inc:29+inc]
#        inc += 6 
    iterator = 1
    dummy = 1
    soup = BeautifulSoup(file)
    list = soup.find_all('table',border='1',limit=1)
    for element in list:
        sublist = element.find_all('tr',recursive=False)
	for subelement in sublist:
            subsublist = subelement.find_all('td',recursive=False)
            for element in subsublist:
                if element.table is None:
		    print iterator,element.get_text()
                    iterator+=1
