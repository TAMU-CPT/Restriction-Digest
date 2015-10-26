#!/usr/bin/env python
"""
DnaDigest digests Bio.Seq objects with enzymes loaded from custom data
structures (usually yaml files)
"""
import re
import yaml
import logging
from pkg_resources import resource_stream
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger()


def __merge_dicts(dict_a, dict_b):
    """Merge dicts where keys map to lists."""
    ret = {}

    for a_dict in (dict_a, dict_b):
        for key in a_dict:
            if key not in ret:
                ret[key] = a_dict[key]
            else:
                ret[key] += a_dict[key]

    return ret


class Enzyme(object):
    """
    A class representing a single enzyme
    """

    DNA_REGEX_TRANSLATIONS = {
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
        self.forward = tuple(forward)  # ('G', 'AATTC')
        self.reverse = tuple(reverse)  # ('CTTAA', 'G')
        assert len(self.forward) == 2
        assert len(self.reverse) == 2
        self.cut_index = len(self.forward[0])

    def _gen_regex_str(self):
        """Generate the expanded regular expression strings for
        forward+backward cuts"""
        reg_f = self.expand_multiple(self.forward[0] + self.forward[1])
        reg_r = self.expand_multiple(self.reverse[0] + self.reverse[1])
        return reg_f, reg_r

    def get_regex(self):
        """Generate the regular expression objects for forward+backward cuts"""
        pre_reg_f, pre_reg_r = self._gen_regex_str()
        reg_f = self.iupac_to_regex(pre_reg_f)
        reg_r = self.iupac_to_regex(pre_reg_r)
        rec_seq_f = re.compile(reg_f, re.IGNORECASE)
        rec_seq_r = re.compile(reg_r, re.IGNORECASE)
        return (rec_seq_f, rec_seq_r)

    def expand_multiple(self, base_str):
        """Sometimes the input sequences contain phrases like N7 which need to
        be expanded to NNNNNNN"""
        match = re.search('(?P<base>[A-Z])(?P<count>[0-9]+)', base_str)
        try:
            # Get position of first match
            base = match.group('base')
            count = int(match.group('count'))

            # Create a fixed string with those bases replaced properly
            replaced = base_str[0:match.start('base')] + \
                base * count + \
                base_str[match.end('count'):]
            # Recurse to replace any more instances of [ACTG][0-9]+
            return self.expand_multiple(replaced)
        except AttributeError:
            return base_str

    def iupac_to_regex(self, recognition_sequence):
        """Replace IUPAC extended DNA alphabet characters with appropriate
        regular experssions from Enzyme.DNA_REGEX_TRANSLATIONS"""
        return ''.join([self.DNA_REGEX_TRANSLATIONS[x] for x in
                        recognition_sequence])

    def __str__(self):
        return "%10s (5' %10s %-10s)" % (self.name, self.forward[0], self.forward[1])

    def digest_iter(self, sequence):
        """Iterator over every cut site in the genome"""
        reg_f, reg_r = self.get_regex()
        LOG.debug('Digest_iter: %s %s', reg_f.pattern, reg_r.pattern)

        for match in reg_f.finditer(str(sequence.seq)):
            cut_location = match.start() + self.cut_index
            LOG.debug('  Match: %s %s %s -> %s', match.start(), match.group(),
                      match.end(), cut_location)
            yield cut_location

        for match in reg_r.finditer(str(sequence.seq)):
            cut_location = match.start() + len(match.group(0)) - self.cut_index
            LOG.debug('  Match: %s %s %s -> %s', match.start(), match.group(),
                      match.end(), cut_location)
            yield cut_location


class DnaDigest(object):
    """Dna Digesting tool with built-in enzyme library."""

    def __init__(self, data_file=None, autoload=True):
        if autoload:
            self.library = self.load_enzyme_data(data_file)
            LOG.debug("Loaded library")
            LOG.debug(self)

    def load_enzyme_data(self, enzyme_data_file):
        """Load enzyme data from a specified path, or from internal rebase.yml
        file is arg[0] is None"""
        if enzyme_data_file is None:
            handle = resource_stream(__name__, 'rebase.yaml')
            return self._load_enzyme_data(handle)
        else:
            with open(enzyme_data_file, 'r') as handle:
                return self._load_enzyme_data(handle)

    @classmethod
    def _load_enzyme_data(cls, data_handle):
        """Actual code to convert rebase.yml file into a usable format
        """
        enzyme_data = {}

        data_structure = yaml.load(data_handle)
        for enzyme_key in data_structure:
            enzyme = data_structure[enzyme_key]
            if len(enzyme['cut'][0]) != len(enzyme['cut'][1]):
                LOG.warning("Cannot use %s; no support for non-matching cuts", enzyme['enzyme'])
            elif len(enzyme['cut']) != 2:
                LOG.warning("Cannot use %s; too many cut sites", enzyme['enzyme'])
            else:
                # Convert
                # d['k'] = ["5' asdfasdf", "3' asdfasdf"]
                # to
                # d['k'] = {"5": "asdfasdf", "3": "asdfasdf" }
                enzyme['cut'] = {
                    x[0]: x[6:-6] for x in enzyme['cut']
                }
            enzyme_data[enzyme_key] = Enzyme(
                name=enzyme_key,
                forward=enzyme['cut']['5'].split('  '),
                reverse=enzyme['cut']['3'].split('  ')
            )
        return enzyme_data

    def digest(self, sequence, enzyme):
        """Digest a signle sequence with a specified enzyme from library"""
        enz = self.library.get(enzyme, None)
        if enz is None:
            raise Exception("Unknown enzyme")

        return self._digest(sequence, enz, did_cut=False, offset=0)

    def multidigest(self, sequences, enzyme):
        """Digest multiple sequences with a specified enzyme from library"""
        enz = self.library.get(enzyme, None)
        if enz is None:
            raise Exception("Unknown enzyme")

        fragments = []
        cut_sites = {}
        did_cut = []

        current_offset = 0
        for sequence in sequences:
            frags, cuts, didc = self._digest(sequence, enz, did_cut=False, offset=current_offset)
            # Update fragment list
            fragments += frags
            cut_sites = __merge_dicts(cut_sites, cuts)
            did_cut.append(didc)
            current_offset += len(sequence)

        return fragments, cut_sites, any(did_cut)

    @classmethod
    def _digest(cls, sequence, enzyme, did_cut=False, offset=0):
        """
        Internal method to digest a sequence with a specified enzyme. Exposed
        publicly in case you wish to use _digest without the enzyme library.
        """
        # TODO: try and make this appx. the length of the cut site, we don't
        # want to have a case where we match TWO times within the wrapped
        # around section
        wrap_around = 15

        cut_sites = []
        fragments = []
        prev_index = 0
        did_cut = did_cut

        if not hasattr(sequence, 'circular'):
            sequence.circular = True

        LOG.debug('Digesting %s with %s. Circ=%s', sequence.id, enzyme, sequence.circular)
        if not sequence.circular:
            # For all the places we hit
            for cut_location in enzyme.digest_iter(sequence):
                LOG.debug("  Cut at %s", cut_location)
                # Store new fragment/cut_site
                fragments.append(sequence[prev_index:cut_location])
                cut_sites.append(cut_location + offset)

                # update our last cut_site
                prev_index = cut_location

            # Append the final fragment
            fragments.append(sequence[prev_index:])
            cut_sites.append(prev_index + offset)
            did_cut = True
        else:  # sequence.circular
            # Add a little bit on the end where it'd "wrap"
            mod_sequence = sequence + sequence[0:wrap_around]

            cut_location = enzyme.digest_iter(mod_sequence).next()
            # If this is the first cut (we'll only ever receive 1+ linear
            # sequences *OR*, **ONE** circular sequence), in order to handle some
            # nasty corner cases more nicely, we'll make the first cut, and
            # then just call/return ourselves with the strand opened at the
            # first cut site as a linear sequence.
            #
            # Linear sequences are handled MUCH more cleanly.
            if cut_location < len(sequence):
                reopened_sequence = sequence[cut_location:] + \
                    sequence[0:cut_location]
                reopened_sequence.circular = False
            else:
                reopened_sequence = mod_sequence[cut_location:] + \
                    mod_sequence[wrap_around:cut_location]
                reopened_sequence.circular = False
            return cls._digest(reopened_sequence, enzyme, did_cut=True, offset=cut_location)

        return fragments, cut_sites, did_cut

    def __str__(self):
        ret = "EnzymeLibrary with %s enzymes\n" % len(self.library)
        for enzyme in self.library:
            ret += "    %s\n" % self.library[enzyme]
        return ret

    def digest_sequence(self, sequence, enzymes):
        """Recursive function to handle digestion with multiple enzymes"""
        enzyme = enzymes[0]
        if len(enzymes) == 1:
            fragments, cut_sites, did_cut = self.digest(sequence, enzyme)
            fixed_cut_sites = {idx: [enzyme] for idx in cut_sites}
            return fragments, fixed_cut_sites, did_cut
        else:
            frags, cuts, didc = self.digest_sequence(self, sequence, enzymes[1:])
            fragments, cut_sites, did_cut = self.multidigest(frags, enzyme)
            return fragments, __merge_dicts(cut_sites, cuts), did_cut or didc
