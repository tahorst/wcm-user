import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from wholecell.analysis.analysis_tools import exportFigure

for setId in ['K', 'L', 'M']:
	reader = csv.reader(open(setId + '.tsv', 'r'), delimiter='\t')
	headers = reader.next()

	divTime = []
	proteomeR = []
	fluxomeR = []

	for line in reader:
		divTime += [float(line[0])]
		proteomeR += [float(line[1])]
		fluxomeR += [float(line[3])]

	fig = plt.figure(figsize = (8.5, 11))
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(divTime, proteomeR, fluxomeR)

	exportFigure(plt, '.', 'shuffle')
	plt.close("all")

	import ipdb; ipdb.set_trace()