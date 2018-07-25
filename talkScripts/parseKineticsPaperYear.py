'''
Script to read tsv of kinetics parameters and determine the year and how many parameters were measured for each paper
Takes input as spreadsheet that needs a pmid and kcat, kM, kI and custom parameter column
Also needs sim_data.cp cPickle file in folder
Outputs file to out.tsv with pmid, year and number of parameters
'''

import urllib
import re
import csv
import os
import json
import time
import cPickle

def getInt(n):
	try:
		return int(n)
	except:
		return n

folder = os.path.dirname(__file__)
inFilename = 'kineticsPapers.tsv'
outFilename = 'out.tsv'

# open tsv with set of kinetics papers and values
kineticsData = csv.reader(open(os.path.join(folder, inFilename), 'r'), delimiter='\t')
headers = kineticsData.next()

pmIdx = headers.index('Pubmed ID')

# list of headers to read values as parameters that have been included in the model
dataHeaders = ['kcat (1/s)', 'kM (uM)', 'kI', 'customParameterConstantValues']
dataIdx = []
for header in dataHeaders:
	dataIdx.append(headers.index(header))

# create dicts for values to be saved
parametersReported = {}
years = {}

# load in years from out file if it exists to prevent excess url calls
if os.path.exists(os.path.join(folder, outFilename)):
	outFile = open(os.path.join(folder, outFilename), 'r')
	outData = csv.reader(outFile, delimiter='\t')
	outData.next()
	for line in outData:
		years[getInt(line[0])] = line[1]
	outFile.close()

print 'Starting at %s' % time.ctime()

for i, line in enumerate(kineticsData):
	if i % 50 == 0:
		print 'row %i' % i
	pmid = getInt(line[pmIdx])

	# determine number of parameters based on length of each data column (each column read as an array)
	nParameters = 0
	for idx in dataIdx:
		if line[idx] != 'null':
			try:
				nParameters += len(json.loads(line[idx]))
			except:
				print "Could not convert to json for %s" % pmid

	parametersReported[pmid] = parametersReported.get(pmid, 0) + nParameters

	if pmid not in years:
		if type(pmid) == str:
			continue
		try:
			# open pubmed site based on pmid of paper
			site = 'https://www.ncbi.nlm.nih.gov/pubmed/%s?report=xml&format=text' % pmid
			url = urllib.urlopen(site)

			# search for the year the article was published
			xml = url.read()
			year = re.search('&lt;DateCreated&gt;.*?&lt;Year&gt;(.*?)&lt;/Year&gt;', xml, re.S).group(1)

			years[pmid] = year
		except KeyboardInterrupt:
			import sys; sys.exit(1)
		except:
			print 'Could not get year for pmid: %s' % pmid

# load sim data to get constraints used in the model
sim_data = cPickle.load(open(os.path.join(folder, "sim_data.cp"), "rb"))

parametersUsed = {}

for constraint in sim_data.process.metabolism.constraintDict.values():
	pmid = constraint['Pubmed ID']

	# 1 parameter for kcat value and additional parameters depending on kM and kI
	nParmaeters = 1 + len(constraint['kM']) + len(constraint['kI'])
	parametersUsed[pmid] = parametersUsed.get(pmid, 0) + nParmaeters

print 'Finished at %s' % time.ctime()

# print out values to new file
output = csv.writer(open(os.path.join(folder, 'out.tsv'), 'w'), delimiter='\t')
output.writerow(['PMID', 'Year Published', 'Parameters Used', 'Parameters Reported'])
for pmid in sorted(years):
	output.writerow([pmid, years[pmid], parametersUsed.get(pmid, 0), parametersReported[pmid]])
import ipdb; ipdb.set_trace()