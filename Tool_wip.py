import pprint
import re
from Bio import SeqIO
import sys
def get_dict():
	enzyme_list = []
	recog_sequence_1 = []
	recog_sequence_2 = []
	cut_list = []
	cut_strand_1=[]
	cut_strand_2 = []
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

def matcher(sequence,enzyme,recognition_sequence):
    mod_seq_string = ''
    for letter in recognition_sequence:
        if letter in 'ATGC':
            mod_seq_string+=letter
        elif letter == 'N':
            mod_seq_string+='.'
        elif letter == 'M':
            mod_seq_string+='[AC]'
        elif letter == 'R':
            mod_seq_string+='[AG]'
        elif letter == 'W':
            mod_seq_string+='[AT]'
        elif letter == 'Y':
            mod_seq_string+='[CT]'
        elif letter == 'S':
            mod_seq_string+='[CG]'
        elif letter == 'K':
            mod_seq_string+='[GT]'
        elif letter == 'H':
            mod_seq_string+='[^G]'
        elif letter == 'B':
            mod_seq_string+='[^A]'
        elif letter == 'V':
            mod_seq_string+='[^T]'
        elif letter == 'D':
            mod_seq_string+='[^C]'
    regex = re.compile(mod_seq_string)
    if len(regex.findall(str(sequence)))!=0:
        return enzyme
    else:
        return 'No'

def string_cutter(sequence,recognition,recog_nucl_index):
    rec_seq = re.compile(recognition)
    match_start = rec_seq.search(str(sequence))
    if match_start is not None:
        return sequence[:rec_seq.search(sequence).start()+recog_nucl_index],sequence[recog_nucl_index+rec_seq.search(sequence).start():]
    else:
        return sequence,'END OF SEQUENCE'

def string_processor(old_fragment_list,recognition,recog_nucl_index):
    new_fragment_list = []
    for fragment in old_fragment_list:
        seq1 = fragment
        seq2 = ''
        while seq2!='END OF SEQUENCE':
            if seq1!=fragment:
                new_fragment_list+=[seq1]
                seq1 = seq2
            seq1,seq2=string_cutter(seq1,recognition,recog_nucl_index)
        if seq2 == 'END OF SEQUENCE':
            new_fragment_list+=[seq1]
    if len(new_fragment_list) in [1,2]  and len(old_fragment_list)==1:
        seq1 = old_fragment_list[0]*2
        store = old_fragment_list[0]
        seq2 = ''
        seq1,seq2 = string_cutter(seq1,recognition,recog_nucl_index)
        single_cleavage_sequence = seq2+seq1[0:(len(store)-len(seq2))]
    return new_fragment_list

def main():
        seqs = [str(record.seq) for record in SeqIO.parse(sys.argv[2],'fasta')]
        print len(seqs)
        can_cleave_list = []
        enzyme_dict = get_dict()
        enzyme_list = str(raw_input('Please enter the names of the restriction enzymes separated only by spaces.')).split(' ')
	template = raw_input('Enter 1 for template strand,0 otherwise.')
        for seq in seqs:
            print 1
            if template ==1:
                seq = seq.reverse_complement()
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
            for enzyme in enzyme_dict:
                if matcher(seq,enzyme,enzyme_dict[enzyme][0][0])!='No':
                    can_cleave_list+= [matcher(seq,enzyme,enzyme_dict[enzyme][0][0])]              
            fragment_list = [seq]
            recognition = [enzyme_dict[enzyme][0][0] for enzyme in enzyme_list]
            recog_nucl_index = [enzyme_dict[enzyme][0][2] for enzyme in enzyme_list]
            for pair in zip(recognition,recog_nucl_index):
                fragment_list = string_processor(fragment_list,pair[0],pair[1])
            print seq,fragment_list
            return fragment_list
            

main()
  
		
