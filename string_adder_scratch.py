string = 'AN'
list = []
#for letter in string:
#    if letter in nucleobase_dict:
#        for value in nucleobase_dict[letter]:
#            list+=[string[0:len(string)-1]+value]
#print list


def string_adder(letter,list):
    new_list = []
    nucleobase_dict = {
        'A':['A'], 'G':['G'], 'C':['C'], 'T':['T'], 'N':['A', 'C', 'G', 'T'],
        'M':['A', 'C'], 'R':['A', 'G'], 'W':['A', 'T'], 'Y':['C', 'T'],
        'S':['C', 'G'], 'K':['G', 'T'], 'H':['A', 'C', 'T'],
        'B':['C', 'G', 'T'], 'V':['A', 'C', 'G'], 'D':['A', 'G', 'T']}
    for element in list:
        for value in nucleobase_dict[letter]:
            new_list+=[element+value]
    return new_list

if __name__ == '__main__':
    string = 'NNN'
    extant_list = ['']
    for letter in string:
        extant_list=string_adder(letter,extant_list)
    if len(extant_list)==64:
        print 'Success'
    else:
        print 'WIP'
