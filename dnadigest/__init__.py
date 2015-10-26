#!/usr/bin/env python
import re
import yaml
import logging
from pkg_resources import resource_stream
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


class Enzyme(object):

    drt = {
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

    def __init__(self, name="", forward=None, reverse=None):
        self.name = name
        self.forward = forward  # ('G', 'AATTC')
        self.reverse = reverse  # ('CTTAA', 'G')
        self.cut_index = len(self.forward[0])

    def generate_regex_str(self, recognition_sequence):
        return ''.join([self.dna_regex_translations[x] for x in
                        recognition_sequence])


class DnaFragment(object):
    LINEAR = 0
    CIRCULAR = 1

    def __init__(self, seqobj, status=LINEAR):
        self.seq = seqobj
        self.status = status


class EnzymeLibrary(object):
    """A library of enzymes, can be used to digest sequences"""

    def __init__(self, data_file=None):
        self.library = self.load_enzyme_data(data_file)
        self.mask = []

    def mask(self, enzymes):
        self.mask = [x for x in enzymes if x in self.library]

    def load_enzyme_data(self, enzyme_data_file):
        """Load enzyme data from a specified path, or from internal rebase.yml
        file is arg[0] is None"""
        if enzyme_data_file is None:
            handle = resource_stream(__name__, 'rebase.yaml')
            return self._load_enzyme_data(handle)
        else:
            with open(enzyme_data_file, 'r') as handle:
                return self._load_enzyme_data(handle)

    def _load_enzyme_data(self, data_handle):
        """Actual code to convert rebase.yml file into a usable format
        """
        ed = {}

        data_structure = yaml.load(data_handle)
        for enzyme_key in data_structure:
            enzyme = data_structure[enzyme_key]
            if len(enzyme['cut'][0]) != len(enzyme['cut'][1]):
                log.warning("Cannot use %s; no support for non-matching cuts" %
                            enzyme['enzyme'])
            elif len(enzyme['cut']) != 2:
                log.warning("Cannot use %s; too many cut sites" %
                            enzyme['enzyme'])
            else:
                # Convert
                # d['k'] = ["5' asdfasdf", "3' asdfasdf"]
                # to
                # d['k'] = {"5": "asdfasdf", "3": "asdfasdf" }
                enzyme['cut'] = {
                    x[0]: x[3:-3] for x in enzyme['cut']
                }
                enzyme['sense_cut_idx'] = self.determine_cut_index(enzyme)

            e = Enzyme(
                name=enzyme_key,
                forward=enzyme['cut']['5'].split('  '),
                reverse=enzyme['cut']['3'].split('  ')
            )
        ed[enzyme_key] = e


class Digest2(object):

    def __init__(self, enzyme_data_file=None):
        """Class to digest DNA strings.

        By default a dataset based on the Wikipedia enzyme list is loaded
        """

        self.enzyme_library = EnzymeLibrary(enzyme_data_file)


class Dnadigest():


    def __merged_iter(cls, *args):
        """Treat multiple iterators as a single iterator
        """
        for arg in args:
            for x in arg:
                yield x

    def string_cutter(self, sequence, recognition_fr, recog_nucl_index,
                      status):
        """Cut a sequence with a 5'+3' cut recognition site
        """
        rec_f_exp = self.expand_multiple(recognition_fr['5'])
        rec_r_exp = self.expand_multiple(recognition_fr['3'])
        rec_seq_f = re.compile(self.generate_regex_str(rec_f_exp), re.IGNORECASE)
        rec_seq_r = re.compile(self.generate_regex_str(rec_r_exp), re.IGNORECASE)

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
                adjusted_recog = \
                    self.__adjust_recog_for_strand(recog_nucl_index, rec_f_exp,
                                                   match.group(0))
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

                    return self.string_cutter(reopened_sequence,
                                              recognition_fr, adjusted_recog,
                                              'linear')

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
                        full_first_fragment = mod_sequence[prev_index:] + \
                            fragments[0]
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
                adjusted_recog = \
                    self.__adjust_recog_for_strand(recog_nucl_index, rec_f_exp,
                                                   match.group(0))
                cut_location = match.start() + adjusted_recog
                fragments.append(sequence[prev_index:cut_location])
                prev_index = cut_location
            fragments.append(sequence[prev_index:])

        # Instead of returning status, if len(fragments) > 1: status='linear'
        return fragments

    def find_cut_sites(self, sequence, enzyme_list, status='circular'):
        """Primarily for use with the drawer() method
        """
        enzymes = self.enzyme_dict_filter(self.enzyme_dict, enzyme_list)
        cut_sites = {}
        for enzyme in enzymes:
            fragments, status, did_cut = self.string_processor(
                [sequence],
                self.enzyme_dict[enzyme]['recognition_sequence'],
                self.enzyme_dict[enzyme]['sense_cut_idx'],
                'circular'
            )
            current_pos = 0
            # We proceed from first to penultimate fragment, marking at site
            # AFTER the current fragment we're examining.
            for fragment in fragments[0:-1]:
                current_pos += len(fragment)

                if current_pos not in cut_sites:
                    cut_sites[current_pos] = []
                cut_sites[current_pos].append(enzyme)
        return cut_sites
    def __adjust_recog_for_strand(self, recog_nucl_index, plus_reference,
                                  matchstr):
        # If the matched group is the plus sense strand, then cut site is FINE
        plus_ref_re = re.compile(self.generate_regex_str(plus_reference))
        if plus_ref_re.match(matchstr):
            return recog_nucl_index
        else:
            # Otherwise, invert it against length of matchstr
            return len(matchstr) - recog_nucl_index

    def string_processor(self, fragment_list, recognition_fr,
                         recog_nucl_index, status):
        new_fragment_list = []
        did_cut = False
        for fragment in fragment_list:
            fragments = self.string_cutter(fragment, recognition_fr,
                                           recog_nucl_index, status)

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

    def process_data(self, seq, cut_with, status='circular'):
        filtered_enzyme_dict = self.enzyme_dict_filter(
            self.enzyme_dict, cut_with)

        fragment_list = [seq]

        cuts = []
        for enzyme in filtered_enzyme_dict:
            log.info('Cutting [%s] with %s' % (','.join(fragment_list),
                                               enzyme))
            log.debug(filtered_enzyme_dict[enzyme])
            (fragment_list, status, did_cut) = \
                self.string_processor(fragment_list,
                                      filtered_enzyme_dict[enzyme]['recognition_sequence'],
                                      filtered_enzyme_dict[enzyme]['sense_cut_idx'],
                                      status)
            if did_cut:
                cuts.append(enzyme)

        return {
            'fragment_list': fragment_list,
            'cut_with': cuts,
            'status': status,
        }
