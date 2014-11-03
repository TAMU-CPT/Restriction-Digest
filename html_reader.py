from bs4 import BeautifulSoup
import sys
import glob
import os

for i in glob.glob(os.path.join(sys.argv[1], '*')):
    with open(i,'r') as file:
        i = 0
        iterator = 1
        j = 0
        enzyme_list = []
        cut_list_1 = []
        cut_list_2 = []
        page = BeautifulSoup(file)
        list_of_tables = page.find_all('table',border='1')
        for table  in list_of_tables:
            list_of_rows  = table.find_all('tr',recursive=False)
            for row  in list_of_rows:
                list_of_cells  = row.find_all('td',recursive=False)
                for cell in list_of_cells:
                    if iterator == 1 and cell.get_text()!= 'Enzyme':
                        enzyme_list += [cell.get_text()]
                    elif iterator == 4 and cell.get_text()!='Recognition sequence':
                        cut_list_1 += [cell.get_text()]
                    elif iterator == 5 and cell.get_text()!='Cut':
                        cut_list_2 += [cell.get_text()]
                    iterator+=1
                iterator = 1
    # print enzyme_list
    '''
    Cleans up the text and writes the list to a file.
    '''
    with open('test.txt','a') as f:
        while i < len(enzyme_list):
            enzyme_string  = enzyme_list[i]
            cut_string_1 =' ['+cut_list_1[i].replace('\n','') + ']'
            cut_string_2 =' ['+cut_list_2[i].replace(u'\xa0','') + ']'
            cut_string_2 = cut_string_2.replace('\n','')
            string = enzyme_string+cut_string_2+'\n'
            #.replace(u'\xa0',u' ')+' '+cut_list_1[i].replace(u'\xa0',u' ').replace(u'\n3',u' ')+' '+cut_list_2[i].replace(u'\xa0',u' ').replace(u'\n3',u' '))
            f.write(string)
            i+=1



