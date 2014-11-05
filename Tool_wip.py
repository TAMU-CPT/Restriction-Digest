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
    if len(rec_seq.findall(str(sequence)))== 0 and status == 'circular':
        working_seq = sequence*2
        start = rec_seq.search(str(sequence))    
        if rec_seq.search(working_seq) is not None and len(rec_seq.findall(working_seq))<=1:
            return working_seq[rec_seq.search(working_seq).start()+recog_nucl_index:len(sequence)]+working_seq[0:rec_seq.search(working_seq).start()+recog_nucl_index],'END OF SEQUENCE','linear'
        else:
            return sequence,'END OF SEQUENCE','circular'
    if len(rec_seq.findall(str(sequence)))== 1 and status == 'circular':
        working_seq = sequence*2
        start = rec_seq.search(str(sequence))
        if rec_seq.search(working_seq) is not None:
            return working_seq[rec_seq.search(working_seq).start()+recog_nucl_index:len(sequence)]+working_seq[0:rec_seq.search(working_seq).start()+recog_nucl_index], working_seq[rec_seq.search(working_seq).start()+recog_nucl_index:len(sequence)]+working_seq[0:rec_seq.search(working_seq).start()+recog_nucl_index],'linear'
        else:
            return sequence,'END OF SEQUENCE','circular'
    if match_start is not None and (len(rec_seq.findall(str(sequence)))!=1 or status!='circular'):
        return sequence[:rec_seq.search(sequence).start()+recog_nucl_index],sequence[recog_nucl_index+rec_seq.search(sequence).start():],'linear'
    else:
        return sequence,'END OF SEQUENCE','linear'

def string_processor(old_fragment_list,recognition,recog_nucl_index,status):
    new_fragment_list = []
    for fragment in old_fragment_list:
        seq1 = fragment
        seq2 = ''
        while seq2!='END OF SEQUENCE':
            if seq1!=fragment and seq2!=seq1:
                new_fragment_list+=[seq1]
                seq1 = seq2
            seq1,seq2,status=string_cutter(seq1,recognition,recog_nucl_index,status)
        if seq2 == 'END OF SEQUENCE':
            new_fragment_list+=[seq1]
    if len(new_fragment_list) in [1,2]  and len(old_fragment_list)==1:
        seq1 = old_fragment_list[0]*2
        store = old_fragment_list[0]
        seq2 = ''
        seq1,seq2,status = string_cutter(seq1,recognition,recog_nucl_index,status)
        single_cleavage_sequence = seq2+seq1[0:(len(store)-len(seq2))]
    return new_fragment_list,status

#def traverse_circular_loop(sequence,recog_length,recog,recog_nucl_index,status):
    #rec_seq = re.compile(recognition)
    #return seq1, 'END OF SEQUENCE',status


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
        for recognition_pair in zip(recognition,recog_nucl_index):
            fragment_list,status = string_processor(fragment_list,recognition_pair[0],recognition_pair[1],status)
        if '' in fragment_list:
            fragment_list.remove('')
        print seq, fragment_list
                 
