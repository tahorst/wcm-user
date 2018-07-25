'''
File to create reaction flat file from ecocyc smart table
Steps to run:
	1. export smart table from ecocyc as spreadsheet (formatted as a tsv, saved as txt)
	2. specify downloaded txt file as FILE - file must be in same directory as script

Note: if exported with frame IDs option then Object ID columns are not needed - will need to adjust idx variables if not the case

TODO:
- add file to run as command line arg
- add export to file with format '%s\t%s\t%s\t%s' % (rid, stoich, reversible, enzymes)
'''

import csv
import re
import os

path = os.path.dirname(__file__)
file = 'ecocycReactions072717.tsv'

reader = csv.reader(open(os.path.join(path, file), 'rb'), delimiter='\t')
headers = reader.next()

# get relavant column ids for information to extract
ridIdx = headers.index('Reactions')
enzymeIdx = headers.index('Enzymes of a reaction')
dirIdx = headers.index('Reaction-Direction')
leftIdx = headers.index('Left')
rightIdx = headers.index('Right')


for line in reader:
	### get reaction id
	rid = line[ridIdx]

	### generate stoich dict
	stoich = {}
	direction = 1
	if 'RIGHT-TO-LEFT' in line[dirIdx]:
		direction = -1

	for left in re.split(' // ', line[leftIdx]):
		# error checking
		if left in stoich:
			print 'Duplicate metabolite for reaction %s' % rid
		# if left.startswith('a '):
		# 	print 'General metabolite in reaction %s' % rid

		stoich[left] = -direction

	for right in re.split(' // ', line[rightIdx]):
		# error checking
		if right in stoich:
			print 'Duplicate metabolite for reaction %s' % rid
		# if right.startswith('a '):
		# 	print 'General metabolite in reaction %s' % rid

		stoich[right] = direction

	# TODO: add check for correct stoich
	# TODO: add compartment info

	### is reversible
	reversible = (line[dirIdx] == 'REVERSIBLE')

	### catalyzed by
	enzymes = re.split(' // ', line[enzymeIdx])
	if len(enzymes) == 1 and len(enzymes[0]) == 0:
		enzymes = []

	import ipdb; ipdb.set_trace()
