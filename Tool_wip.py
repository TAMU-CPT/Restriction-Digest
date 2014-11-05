import pprint
import re
from Bio import SeqIO
import sys
import yaml
import argparse


def get_dict():
    with open('enzyme_data.yaml') as handle:
        data_structure = yaml.load(handle)

    tmp_corrected = {}
    for enzyme in data_structure:
        #'ZhoI': (['ATCGAT', '--AT  CGAT---', ''], ['TAGCTA', '---TAGC
        # TA---', '']),
        try:
            tmp_corrected[enzyme['enzyme']] = (
                [
                    enzyme['recognition_sequence'][0].split()[1],
                    enzyme['cut'][0].split()[1],
                    ''
                ],
                [
                    enzyme['recognition_sequence'][1].split()[1],
                    enzyme['cut'][1].split()[1],
                    ''
                ],
            )
        except Exception, e:
            # These enzymes will need to be corrected, possibly on wikipedia.
            #print e
            #print enzyme
            #Dummy line to prevent the error message from popping up every time I run the program
            x=2
    return tmp_corrected

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

def string_cutter(sequence,recognition,recog_nucl_index,status):
    rec_seq = re.compile(recognition)
    match_start = rec_seq.search(str(sequence))
    if len(rec_seq.findall(str(sequence))) == 1 and status == 'circular':
        return traverse_circular_loop(sequence,len(recognition),recognition,status)
    if match_start is not None:
        return sequence[:rec_seq.search(sequence).start()+recog_nucl_index],sequence[recog_nucl_index+rec_seq.search(sequence).start():]
    else:
        return sequence,'END OF SEQUENCE'

def string_processor(old_fragment_list,recognition,recog_nucl_index,status):
    new_fragment_list = []
    for fragment in old_fragment_list:
        seq1 = fragment
        seq2 = ''
        while seq2!='END OF SEQUENCE':
            if seq1!=fragment:
                new_fragment_list+=[seq1]
                seq1 = seq2
            print 'Sequence:',seq1
            print 'Recognition:',recognition
            print 'Index:',recog_nucl_index
            print 'Status:',status
            #seq1,seq2,status=string_cutter(seq1,recognition,recog_nucl_index,status)
            print string_cutter(seq1,recognition,recog_nucl_index,status)
        if seq2 == 'END OF SEQUENCE':
            new_fragment_list+=[seq1]
    if len(new_fragment_list) in [1,2]  and len(old_fragment_list)==1:
        seq1 = old_fragment_list[0]*2
        store = old_fragment_list[0]
        seq2 = ''
        seq1,seq2 = string_cutter(seq1,recognition,recog_nucl_index)
        single_cleavage_sequence = seq2+seq1[0:(len(store)-len(seq2))]
    return new_fragment_list

def traverse_circular_loop(string,step,recog,status):
    working_str = string*2
    cut = [0]
    start = 0
    end = step
    string2 = '  '
    while string2[1]!=working_str[0]:
        string2 = working_str[start:end]
        start+=1
        if string2 == recog:
            cut[0] = end
        end+=1
    print string
    if cut[0]>= len(string):
        cut[0] = cut[0]%len(string)
    if cut[0]!=0:
        status = 'linear'
        return string[cut:]+string[0:cut],' ',status
    else:
        return string,' ',status

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restriction Digest Tool')
    # Input filename, required
    parser.add_argument('file', help='Input fasta genome(s)')
    # A single string with default value of 'enzyme_data.yaml'
    parser.add_argument('--data', help='Enzyme cut site dataset', default='enzyme_data.yaml')
    # A list of one or more strings, at the end
    parser.add_argument('enzyme', metavar='E', type=str, nargs='+', help='Space separated list of enzymes')
    args = parser.parse_args()

    seqs = [str(record.seq) for record in SeqIO.parse(args.file,'fasta')]
    can_cleave_list = []
    enzyme_dict = get_dict()
    # For now, just write against plus strand, this is a minor issue that can
    # be corrected later: This program should function against BOTH strands
    # without asking user, as that's the biological reality. If only we could
    # selectively cut against only one strand...

    #template = raw_input('Enter 1 for template strand,0 otherwise.')
    for seq in seqs:
        status = 'circular'
        #if template ==1:
            #seq = seq.reverse_complement()
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
        recognition = [enzyme_dict[enzyme][0][0] for enzyme in args.enzyme]
        recog_nucl_index = [enzyme_dict[enzyme][0][2] for enzyme in args.enzyme]
        for pair in zip(recognition,recog_nucl_index):
            fragment_list = string_processor(fragment_list,pair[0],pair[1],status)
        print seq, fragment_list
         
