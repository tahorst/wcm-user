#! /usr/bin/env python

"""
Template file for new python files.
"""

import argparse
import os
import time


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))


def main() -> None:
	"""
	Main function
	"""

	return

def parse_args() -> argparse.Namespace:
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('-a', '--abc',
		type=int,
		default=1,
		help='help')
	parser.add_argument('--true',
		action='store_true',
		help='help')

	return parser.parse_args()


if __name__ == '__main__':
	start = time.time()

	args = parse_args()

	main()

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
