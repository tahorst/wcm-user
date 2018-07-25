import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from memory_profiler import profile
import os
import gc

parent = os.path.dirname(__file__)

@profile
def plot(idx):
	plt.figure()

	plt.plot(range(1000), range(1000))
	plt.savefig(os.path.join(parent, 'png', '%s.png' % idx))
	plt.close('all')


@profile
def plot_oo(idx):
	fig, ax = plt.subplots()

	ax.plot(range(1000), range(1000))
	plt.savefig(os.path.join(parent, 'png', '%s.png' % idx))
	plt.close('all')


@profile
def plot_gc(idx):
	plt.figure()

	plt.plot(range(1000), range(1000))
	plt.savefig(os.path.join(parent, 'png', '%s.png' % idx))
	plt.close('all')

	gc.collect()


if __name__ == '__main__':
	for i in range(2000):
		plot(i)

	for i in range(3):
		plot_oo(i)

	for i in range(3):
		plot_gc(i)