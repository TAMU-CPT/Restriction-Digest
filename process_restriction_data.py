import pprint
tracker = 0
enzyme_list = []
recog_sequence_1 = []
recog_sequence_2 = []
cut_list = []
cut_strand_1=[]
cut_strand_2 = []
cut_start_list = []
info_dict = {}
m=0
with open('test.txt','r') as file:
    for line in file.readlines():
	i=0
        while line[i]!=' ':
		i+=1
        j=-1
        while line[j]!='[':
                j-=1
        cut_list += [line[j:]]
	enzyme_list += [line[0:i]]
    for cut in cut_list:
        k=-1
        l=0
        while cut[k]!='3':
            k-=1
        while cut[l]!='3':
            l+=1
        cut_strand_1+=[cut[5:l]]
        cut_strand_2+=[cut[k+2:len(cut)-4]]
    for cut in cut_strand_1:
        string= ''
        for letter in cut:
            if letter!=' ' and letter!='-':
                string+=letter
        recog_sequence_1+=[string]
    for cut in cut_strand_2:
        string= ''
        for letter in cut:
            if letter!=' ' and letter!='-':
                string+=letter
        recog_sequence_2+=[string]
    for enzyme in enzyme_list:
        info_dict[enzyme] = ((recog_sequence_1[m],cut_strand_1[m]),(recog_sequence_2[m],cut_strand_2[m]))
        m+=1
    pprint.pprint(info_dict)
    
        
		
