import re
import yaml
import logging
logging.basicConfig()
log = logging.getLogger()


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

        # These provide translations between the ambiguity codes.
        self.complement_regex = {
            'A': 'T',
            'C': 'G',
            'N': 'N',
            'H': 'D',
            'B': 'V',
            'M': 'K',
            'R': 'Y',
            'W': 'S',
        }
        # Less typing
        for k in self.complement_regex.keys():
            self.complement_regex[self.complement_regex[k]] = k

    def get_dict(self, data_file):
        with open(data_file) as handle:
            data_structure = yaml.load(handle)

        tmp_corrected = {}
        for enzyme in data_structure:
            if len(enzyme['cut'][0]) != len(enzyme['cut'][1]):
                log.warning("Cannot use %s; no support for non-matching cuts" % enzyme['enzyme'])
            elif len(enzyme['cut']) != 2:
                log.warning("Cannot use %s; too many cut sites" % enzyme['enzyme'])
            else:

                # Convert
                # d['k'] = ["5' asdfasdf", "3' asdfasdf"]
                # to
                # d['k'] = {"5": "asdfasdf", "3": "asdfasdf" }
                enzyme['recognition_sequence'] = {
                    x[0]: x[3:] for x in enzyme['recognition_sequence']
                }
                enzyme['cut'] = {
                    x[0]: x[3:-3] for x in enzyme['cut']
                }
                enzyme['sense_cut_idx'] = self.determine_cut_index(enzyme)


                tmp_corrected[enzyme['enzyme']] = enzyme
        return tmp_corrected

    @classmethod
    def determine_cut_index(cls, enzyme):
        return enzyme['cut']['5'].strip('-').index(' ')

    def generate_regex_str(self, recognition_sequence):
        return ''.join([self.dna_regex_translations[x] for x in
                        recognition_sequence])

    def matcher(self, sequence, recognition_sequence):
        regex = re.compile(self.generate_regex_str(recognition_sequence))
        return len(regex.findall(sequence)) != 0

    def __merged_iter(cls, *args):
        """Treat multiple iterators as a single iterator
        """
        for arg in args:
            for x in arg:
                yield x

    def string_cutter(self, sequence, recognition_fr, recog_nucl_index, status):
        """Cut a sequence with a 5'+3' cut recognition site
        """
        rec_f_exp = self.expand_multiple(recognition_fr['5'])
        rec_r_exp = self.expand_multiple(recognition_fr['3'])
        rec_seq_f = re.compile(self.generate_regex_str(rec_f_exp))
        rec_seq_r = re.compile(self.generate_regex_str(rec_r_exp))

        # TODO: try and make this appx. the length of the cut site, we don't
        # want to have a case where we match TWO times within the wrapped
        # around section
        wrap_around = 15

        fragments = []
        prev_index = 0

        if status == 'circular':
            # Add a little bit on the end where it'd "wrap"
            mod_sequence = sequence + sequence[0:wrap_around]
            match_list = self.__merged_iter(
                rec_seq_f.finditer(mod_sequence),
                rec_seq_r.finditer(mod_sequence)
            )
            # Track where our first cut was made
            first_cut = None
            # Cleanup for some corner cases
            remove_first_fragment = False
            for match in match_list:
                adjusted_recog = self.__adjust_recog_for_strand(recog_nucl_index, rec_f_exp, match.group(0))
                cut_location = match.start() + adjusted_recog
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

                    return self.string_cutter(reopened_sequence, recognition_fr,
                                              adjusted_recog, 'linear')

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
            match_list = self.__merged_iter(
                rec_seq_f.finditer(sequence),
                rec_seq_r.finditer(sequence)
            )
            for match in match_list:
                adjusted_recog = self.__adjust_recog_for_strand(recog_nucl_index, rec_f_exp, match.group(0))
                cut_location = match.start() + adjusted_recog
                fragments.append(sequence[prev_index:cut_location])
                prev_index = cut_location
            fragments.append(sequence[prev_index:])

        # Instead of returning status, if len(fragments) > 1: status='linear'
        return fragments

    def __adjust_recog_for_strand(self, recog_nucl_index, plus_reference, matchstr):
        #log.debug("RNI: %s, PR: %s, M: %s" % (recog_nucl_index, plus_reference, matchstr))
        # If the matched group is the plus sense strand, then cut site is FINE
        plus_ref_re = re.compile(self.generate_regex_str(plus_reference))
        #log.debug('%s %s %s' % (plus_reference, plus_ref_re.match(matchstr), matchstr))
        if plus_ref_re.match(matchstr):
            return recog_nucl_index
        else:
            # Otherwise, invert it against length of matchstr
            return len(matchstr) - recog_nucl_index

    def string_processor(self, fragment_list, recognition_fr, recog_nucl_index, status):
        new_fragment_list = []
        did_cut = False
        log.debug("SP: FL: %s, RF: %s, RNI: %s, ST: %s" % (fragment_list, recognition_fr, recog_nucl_index, status))
        for fragment in fragment_list:
            fragments = self.string_cutter(fragment, recognition_fr, recog_nucl_index, status)
            log.info("Fragments: [%s]" % ','.join(fragments))

            if status == 'circular' and len(fragments) > 0:
                status = 'linear'

            if len(fragments) > 0:
                did_cut = True

            new_fragment_list += fragments

        # Ensure we return a complete fragment list and not empty
        if len(new_fragment_list) == 0:
            return fragment_list, status, False
        else:
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

    def enzyme_dict_filter(self, data, cut_list):
        # TODO: need to include isoscizomers, but current data structure
        # doesn't allow for that.
        #
        # For the time being, just remove all enzymes that the user didn't
        # request
        good = {}
        for enzyme in data:
            if enzyme in cut_list:
                good[enzyme] = data[enzyme]
            elif 'isoscizomers' in data[enzyme]:
                for iso in data[enzyme]['isoscizomers']:
                    if iso in cut_list:
                        good[enzyme] = data[enzyme]
                        continue
        return good

    def process_data(self, seq, enzyme_dict, cut_with, status='circular'):
        # For now, just write against plus strand, this is a minor issue that can
        # be corrected later: This program should function against BOTH strands
        # without asking user, as that's the biological reality. If only we could
        # selectively cut against only one strand...

        # Filter for requested enzymes so we aren't searching and processing
        # HUGE amounts of data
        enzyme_dict = self.enzyme_dict_filter(enzyme_dict, cut_with)

        fragment_list = [seq]

        for enzyme in enzyme_dict:
            #print 'Cutting ', fragment_list, 'with',  enzyme
            log.info('Cutting [%s] with %s' % (','.join(fragment_list), enzyme))
            (fragment_list, status, did_cut) = \
                self.string_processor(fragment_list,
                                        enzyme_dict[enzyme]['recognition_sequence'],
                                        enzyme_dict[enzyme]['sense_cut_idx'],
                                        status)
        return fragment_list, status
