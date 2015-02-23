import re
import yaml


class Dnadigest():
    def __init__(self):
        self.data = ''
        self.dna_regex_translations = {
            'A': 'A',
            'T': 'T',
            'C': 'C',
            'G': 'G',
            'N': '.',
            'M': '[AC]',
            'R': '[AG]',
            'W': '[AT]',
            'Y': '[CT]',
            'S': '[CG]',
            'K': '[GT]',
            'H': '[^G]',
            'B': '[^A]',
            'V': '[^T]',
            'D': '[^C]',
        }

    def get_dict(self, data_file):
        with open(data_file) as handle:
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
            except:
                # These enzymes will need to be corrected, possibly on wikipedia.
                pass

        # GLobal data which requires global processing MUST be done on load,
        # not for each processing loop
        for enzyme in tmp_corrected:
            tmp_corrected[enzyme] = self.__find_cut_site(tmp_corrected[enzyme])

        return tmp_corrected

    def generate_regex_str(self, recognition_sequence):
        return ''.join([self.dna_regex_translations[x] for x in recognition_sequence])

    def matcher(self, sequence, recognition_sequence):
        regex = re.compile(self.generate_regex_str(recognition_sequence))
        return  len(regex.findall(sequence)) != 0

    def string_cutter(self, sequence, recognition, recog_nucl_index, status):
        """Attempt to cut a piece of DNA with a specific enzyme
        """
        # TODO: Refactor!!!!
        # This code "smells" really bad.
        rec_seq = re.compile(self.generate_regex_str(recognition))

        if status == 'circular':
            # Add a little bit on the end where it'd "wrap"
            match_list = rec_seq.search(sequence + sequence[0:100])
        else:
            match_list = rec_seq.search(sequence)

        print match_list


        if len(rec_seq.findall(str(sequence)))== 0 and status == 'circular':
            working_seq = sequence*2
            if rec_seq.search(working_seq) is not None and len(rec_seq.findall(working_seq))<=1:
                return working_seq[rec_seq.search(working_seq).start()+recog_nucl_index:len(sequence)]+working_seq[0:rec_seq.search(working_seq).start()+recog_nucl_index],'END OF SEQUENCE','linear',rec_seq.search(working_seq).start()+recog_nucl_index
            else:
                return sequence,'END OF SEQUENCE','circular',''
        if len(rec_seq.findall(str(sequence)))== 1 and status == 'circular':
            working_seq = sequence*2
            if rec_seq.search(working_seq) is not None:
                return working_seq[rec_seq.search(working_seq).start()+recog_nucl_index:len(sequence)]+working_seq[0:rec_seq.search(working_seq).start()+recog_nucl_index], working_seq[rec_seq.search(working_seq).start()+recog_nucl_index:len(sequence)]+working_seq[0:rec_seq.search(working_seq).start()+recog_nucl_index],'linear',rec_seq.search(working_seq).start()+recog_nucl_index
            else:
                return sequence,'END OF SEQUENCE','circular',''
        if match_start is not None and (len(rec_seq.findall(str(sequence)))!=1 or status!='circular'):
            return (
                sequence[:rec_seq.search(sequence).start()+recog_nucl_index],
                sequence[recog_nucl_index+rec_seq.search(sequence).start():],
                'linear',
                rec_seq.search(sequence).start()+recog_nucl_index
                )
        else:
            return sequence,'END OF SEQUENCE','linear',''

    def string_processor(self, old_fragment_list, recognition, recog_nucl_index, status):
        new_fragment_list = []
        line_marker_list = []
        line_marker_init = 0

        for fragment in old_fragment_list:
            print "SP: fragment: " + fragment
            print "Cutting: %s %s %s" % (recognition, recog_nucl_index, status)
            print self.string_cutter(fragment, recognition, recog_nucl_index, status)
            seq1 = fragment
            seq2 = ''
            while seq2!='END OF SEQUENCE':
                if seq1!=fragment and seq2!=seq1:
                    new_fragment_list+=[seq1]
                    line_marker_list+=[int(line_marker)+int(line_marker_init)]
                    line_marker_init+=int(line_marker)
                    seq1 = seq2
                seq1,seq2,status,line_marker=self.string_cutter(seq1,recognition,recog_nucl_index,status)
            if seq2 == 'END OF SEQUENCE':
                new_fragment_list+=[seq1]
        if len(new_fragment_list) in [1,2]  and len(old_fragment_list)==1:
            seq1 = old_fragment_list[0]*2
            store = old_fragment_list[0]
            seq2 = ''
            seq1,seq2,status,line_marker = self.string_cutter(seq1,recognition,recog_nucl_index,status)
            single_cleavage_sequence = seq2+seq1[0:(len(store)-len(seq2))]
        return new_fragment_list,status,line_marker_list

    #def traverse_circular_loop(sequence,recog_length,recog,recog_nucl_index,status):
        #rec_seq = re.compile(recognition)
        #return seq1, 'END OF SEQUENCE',status


    def expand_multiple(self, base_str):
        m = re.search('(?P<base>[A-Z])(?P<count>[0-9]+)', base_str)
        try:
            # Get position of first match
            base = m.group('base')
            count = int(m.group('count'))

            # Create a fixed string with those bases replaced properly
            replaced = base_str[0:m.start('base')] + \
                base * count + \
                base_str[m.end('count'):]
            # Recurse to replace any more instances of [ACTG][0-9]+
            return self.expand_multiple(replaced)
        except AttributeError:
            return base_str

    def __find_cut_site(self, enzyme_dict):
        for recogsite, datalist in enumerate(enzyme_dict):
            for i, element in enumerate(enzyme_dict[recogsite]):
                # Expand stuff like N6
                enzyme_dict[recogsite][i] = \
                    self.expand_multiple(element)
        # TODO, check that data doesn't include any ------ACTG-, but it
        # SHOULDN'T, and that there aren't spaces in it...
        enzyme_dict[recogsite][2] = enzyme_dict[recogsite][1].strip().count('-')
        return enzyme_dict[recogsite]

    def enzyme_dict_filter(self, data, cut_list):
        # TODO: need to include isoscizomers, but current data structure
        # doesn't allow for that.
        #
        # For the time being, just remove all enzymes that the user didn't
        # request
        return {x: data[x] for x in data if x in cut_list}

    def process_data(self, seqs, enzyme_dict, cut_with):
        status = 'circular'
        # For now, just write against plus strand, this is a minor issue that can
        # be corrected later: This program should function against BOTH strands
        # without asking user, as that's the biological reality. If only we could
        # selectively cut against only one strand...

        # Filter for requested enzymes so we aren't searching and processing
        # HUGE amounts of data
        enzyme_dict = self.enzyme_dict_filter(enzyme_dict, cut_with)

        assoc_enzyme_list = []
        for seq in seqs:
            fragment_list = [seq]

            for enzyme in enzyme_dict:
                (fragment_list, status, line_marker_list) = self.string_processor(fragment_list,
                                                                                  enzyme_dict[enzyme][0],
                                                                                  enzyme_dict[enzyme][2],
                                                                                  status)
                print "FL: ", fragment_list
                print "Status: ", status
                print "LML: ", line_marker_list

                assoc_enzyme_list.append([cut_with[q] for mark in line_marker_list])

        return fragment_list,assoc_enzyme_list,line_marker_list,len(seq)
