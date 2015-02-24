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
            # 'ZhoI': (['ATCGAT', '--AT  CGAT---', ''], ['TAGCTA', '---TAGC
            #  TA---', '']),
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
        return len(regex.findall(sequence)) != 0

    def string_cutter(self, sequence, recognition, recog_nucl_index, status):
        """Attempt to cut a piece of DNA with a specific enzyme
        """
        rec_seq = re.compile(self.generate_regex_str(recognition))

        # TODO: try and make this appx. the length of the cut site, we don't
        # want to have a case where we match TWO times within the wrapped
        # around section
        wrap_around = 15

        fragments = []
        prev_index = 0

        if status == 'circular':
            # Add a little bit on the end where it'd "wrap"
            mod_sequence = sequence + sequence[0:wrap_around]
            match_list = rec_seq.finditer(mod_sequence)
            # Track where our first cut was made
            first_cut = None
            # Cleanup for some corner cases
            remove_first_fragment = False
            for match in match_list:
                cut_location = match.start() + recog_nucl_index
                if first_cut is None:
                    # If this is the first cut, in order to handle some nasty
                    # corner cases more nicely, just recursively call ourselves
                    # with the strand opened at the first cut site.
                    if cut_location < len(sequence):
                        reopened_sequence = sequence[cut_location:] + \
                            sequence[0:cut_location]
                    else:
                        reopened_sequence = mod_sequence[cut_location:] + \
                            mod_sequence[wrap_around:cut_location]

                    return self.string_cutter(reopened_sequence, recognition,
                                              recog_nucl_index, 'linear')

                # If this is a "normal" cut, append the new fragment from the
                # previous cut site to here
                remove_first_fragment = True
                if cut_location < len(sequence):
                    fragments.append(mod_sequence[prev_index:cut_location])
                    prev_index = cut_location
                else:
                    # This is not a normal cut, i.e. in the wrapped sequenc
                    # This case is a bit complicated:
                    # - need to add the correct fragment
                    # - need to removeleading characters from the first
                    #   fragment (and ensure it wasn't detected there too)
                    if first_cut == cut_location - len(sequence):
                        # First cut was in the same position as this, so we
                        # delete first fragment and trim up to this cut
                        # location here.
                        fragments.append(mod_sequence[prev_index:cut_location])
                        break
                    else:
                        # This cut was NOT caught by the first cut, so this
                        # means that there's some serious overlap and we cannot
                        # delete the first fragment.
                        #
                        # This is a REALLY unpleasant case.

                        # Get the full first fragment by taking the first
                        # fragment with the "latest" sequence, not including
                        # the wrap around
                        full_first_fragment = mod_sequence[prev_index:] + fragments[0]
                        remapped_cut_location = cut_location - prev_index
                        fragments.append(full_first_fragment[0:remapped_cut_location])
                        fragments.append(full_first_fragment[remapped_cut_location:])

            if remove_first_fragment and len(fragments) > 1:
                del fragments[0]
        else:
            match_list = rec_seq.finditer(sequence)
            for match in match_list:
                cut_location = match.start() + recog_nucl_index
                fragments.append(sequence[prev_index:cut_location])
                prev_index = cut_location
            fragments.append(sequence[prev_index:])

        # Instead of returning status, if len(fragments) > 1: status='linear'
        return fragments

    def string_processor(self, old_fragment_list, recognition, recog_nucl_index, status):
        new_fragment_list = []

        did_cut = False

        for fragment in old_fragment_list:
            print "SP: fragment: " + fragment
            print "Cutting: %s %s %s" % (recognition, recog_nucl_index, status)
            fragments = self.string_cutter(fragment, recognition, recog_nucl_index, status)
            print fragments

            if status == 'circular' and len(fragments) > 0:
                status = 'linear'

            if len(fragments) > 0:
                did_cut = True

            new_fragment_list += fragments

        return new_fragment_list, status, did_cut

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

        for seq in seqs:
            fragment_list = [seq]

            for enzyme in enzyme_dict:
                (fragment_list, status, did_cut) = \
                    self.string_processor(fragment_list,
                                          enzyme_dict[enzyme][0],
                                          enzyme_dict[enzyme][2],
                                          status)

                print "FL: ", fragment_list
                print "Status: ", status
                print "DC: ", did_cut

        try:
            del fragment_list[0]
        except:
            pass

        return fragment_list
