import pprint
def get_dict():
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
	with open('Enzyme_Data.txt','r') as file:
	    for line in file.readlines():
		i=0
	        while line[i]!=' ':
			i+=1
	        j=-1
	        while line[j]!='[':
	                j-=1
	        cut_list += [line[j:]]
		enzyme_list += [line[0:i]]
	   # print 'No. of Enzymes:', len(enzyme_list)
	   # print 'No. of Cuts:',cut_list
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
	        info_dict[enzyme] = ([recog_sequence_1[m].strip(),cut_strand_1[m].strip(),''],[recog_sequence_2[m].strip(),cut_strand_2[m].strip(),''])
	        m+=1
	    return info_dict

def get_seq(file):
	'''
	Takes a FASTA file and returns the sequence to be cleaved.
	'''
	return 'ATGATGTTGTGTGTCGCGCCGCGAATGTGTTGACACTGTACATCGATCAGCTAGCTAGCTAGCTAACATGCTGTGTAATATTATATGCGATGCATGC'

def main():
        enzymes_that_can_cleave_before_edit = 172
        i=0
	import sys
        enzyme_dict = get_dict()
        seq = get_seq(open(sys.argv[1],'r'))
        enzyme_list = str(raw_input('Please enter the names of the restriction enzymes separated only by spaces.')).split(' ')
	#Code will be revised to include template strand if necessary.
        for enzyme in enzyme_dict:
            for list in range(2):
                for element in range(2):
                    num_list = ['0','1','2','3','4','5','6','7','8','9']
                    num_string = '0'
                    new_seq = ''
                    for letter in enzyme_dict[enzyme][list][element]:
                        if letter in num_list:
                            num_string+=letter
                        else:
                            if int(num_string)!=0:
                                new_seq+=new_seq[-1]*(int(num_string)-1)+letter
                                num_string = '0'
                            else:
                                new_seq+=letter
                    enzyme_dict[enzyme][list][element]=new_seq
                cut_pos = 0
                for letter in enzyme_dict[enzyme][list][1]:
                    if letter != ' ' and letter!= '-':
                        cut_pos+=1
                    if letter == ' ':
                        break
                enzyme_dict[enzyme][list][2]=cut_pos
        pprint.pprint(enzyme_dict)                   
        nucleobase_dict = {'N':['A','C','G','T'],'M':['A','C'],'R':['A','G'],'W':['A','T'],'Y':['C','T'],'S':['C','G'],'K':['G','T'],'H':['A','C','T'],'B':['C','G','T'],'V':['A','C','G'],'D':['A','G','T']}

main()
        
		
