from bs4 import BeautifulSoup
import os
import yaml
import string
import urllib2

def f(s):
    # Remove non-printable characters/bad characters, and encode down to ascii.
    # Probably only /really/ need to do second step, but just in case.
    return ''.join(c for c in s if c in string.printable).encode('ascii', 'ignore')

enzyme_list = []

#TODO: make a request to the base page and get the list from that page
#for i in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
          #'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'Ba-Bc', 'Bd-Bp',
          #'Bsa-Bso', 'Bsp', 'Bsr-Bss', 'Bst', 'Bsu-Bv']:
for i in ['A', 'C-D', 'E-F', 'G-K', 'L-N', 'O-R', 'S', 'T-Z',
          'Ba-Bc', 'Bd-Bp', 'Bsa-Bso', 'Bsp-Bss', 'Bst-Bv']:

    html_file = "%s.html" % i

    if not os.path.exists(html_file):
        try:
            response = urllib2.urlopen(
                "http://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_%s" % i)
            html_page = response.read()
            # Store the page to file
            with open(html_file, 'w') as handle:
                handle.write(html_page)
        except:
            print "Error accessing %s" % i

    page = BeautifulSoup(open(html_file))
    list_of_tables = page.find_all('table',border='1')
    for table  in list_of_tables:
        list_of_rows  = table.find_all('tr',recursive=False)
        # Skip first line, header
        for row  in list_of_rows[1:]:
            list_of_cells = row.find_all('td',recursive=False)
            # Cells are either newline separated, or comma separated. We
            # need to fix the string data, then split on '\n', as we would
            # rather have recognition_sequence and cut_site be two objects
            # in a list, rather than a single string.
            enzyme_list.append({
                'enzyme': f(list_of_cells[0].get_text()),
                'recognition_sequence': f(list_of_cells[3].get_text()).split('\n'),
                'cut': f(list_of_cells[4].get_text()).split('\n'),
                # These are comma separated, and surrounded by whitespace,
                # thanks wikipeda. So, run .strip() on all to clean them up
                'isoscizomers': [e.strip() for e in f(list_of_cells[5].get_text()).split(',')],
            })

with open('enzyme_data.yaml', 'w') as handle:
    handle.write(yaml.dump(enzyme_list, default_flow_style=False))
