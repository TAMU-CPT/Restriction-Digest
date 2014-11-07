from bs4 import BeautifulSoup
import sys
import glob
import os
import yaml
import string

def f(s):
    # Remove non-printable characters/bad characters, and encode down to ascii.
    # Probably only /really/ need to do second step, but just in case.
    return ''.join(c for c in s if c in string.printable).encode('ascii', 'ignore')

enzyme_list = []
for i in glob.glob(os.path.join(sys.argv[1], '*')):
    print i
    with open(i,'r') as file:
        page = BeautifulSoup(file)
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
