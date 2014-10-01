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
#    iterator = 1
#    dummy = 1
#    soup = BeautifulSoup(file)
#    list = soup.find_all('table',border='1',limit=1)
#    for element in list:
#        sublist = element.find_all('tr',recursive=False)
#	for subelement in sublist:
#            subsublist = subelement.find_all('td',recursive=False)
#            for element in subsublist:
#                if element.table is None:
#		    print iterator,element.get_text()
#                    iterator+=1
#    iterator = 1
#    dummy = 1
#    soup = BeautifulSoup(file)
#    list = soup.find_all('table',border='1',limit=1)
#    for element in list:
#        sublist = element.find_all('tr',recursive=False)
#        for subelement in sublist:
#            subsublist = subelement.find_all('td',recursive=False)
#            for element in subsublist:
#                if element.table is None:
#                    print iterator,element.get_text()
#                    iterator+=1
#                else:
#                    print iterator
#                    iterator+=1



    iterator = 1
    dummy = 1
    enzyme_list = []
    cut_list_1 = []
    cut_list_2 = []
    soup = BeautifulSoup(file)
    list = soup.find_all('table',border='1',limit=1)
    for element in list:
        sublist = element.find_all('tr',recursive=False)
        for subelement in sublist:
            subsublist = subelement.find_all('td',recursive=False)
            for element in subsublist:
                if iterator%6==1:
                    enzyme_list+=[element.get_text()]
                elif iterator%6==4:
                    cut_list_1+=[element.get_text()]
                elif iterator%6==5:
                    cut_list_2+=[element.get_text()]
                iterator+=1
    print enzyme_list
    print cut_list_1
    print cut_list_2
#list = soup.find_all('table',border='1',limit=1)
#for element in list:
#    sublist = element.find_all('tr',recursive=False)
#    for subelement in sublist:
#        subsublist = subelement.find_all('td',recursive=False)
#print subsublist[0:20]
