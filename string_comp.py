import re
import sys
def sys_type_checker():
    print sys.argv[1]
    print type(sys.argv[1])

def string_cutter(sequence,recognition,recog_nucl_index):
    rec_seq = re.compile(recognition)
    match_start = rec_seq.search(sequence)
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
        return [single_cleavage_sequence]
    return new_fragment_list


def main():
    fragment_list = ['ABCDEFG']
    recognition = ['GAB']
    recog_nucl_index = [3]
    #seq1 = sequence
    #seq2 = ''
    for pair in zip(recognition,recog_nucl_index):
        fragment_list = string_processor(fragment_list,pair[0],pair[1])
        print fragment_list

main()

