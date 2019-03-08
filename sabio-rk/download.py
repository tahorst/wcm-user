'''
Script to get all kinetic information related to E coli from SABIO-RK.
Script adapted from example at http://sabio.h-its.org/layouts/content/docuRESTfulWeb/searchPython.gsp

Saves data to ecoli-kinetics.tsv in the script directory
'''

import requests
import os

ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'
PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ecoli-kinetics.tsv')

# ask SABIO-RK for all EntryIDs matching a query
query = {'format':'txt', 'q':'Organism:"Escherichia coli" AND EnzymeType:"wildtype"'}

# make GET request
request = requests.get(ENTRYID_QUERY_URL, params = query)
request.raise_for_status() # raise if 404 error

# each entry is reported on a new line
entryIDs = [int(x) for x in request.text.strip().split('\n')]
print('%d matching entries found.' % len(entryIDs))

# encode next request, for parameter data given entry IDs
data_field = {'entryIDs[]': entryIDs}
query = {'format':'tsv', 'fields[]':['EntryID', 'PubMedID', 'Organism', 'UniprotID', 'Enzymename', 'ECNumber', 'Parameter', 'ReactionEquation']}

# make POST request
request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
request.raise_for_status()

# results to tsv file
with open(OUTPUT_FILE, 'w') as f:
	f.write(request.text)
