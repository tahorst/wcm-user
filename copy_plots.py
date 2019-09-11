#! /usr/bin/env python
'''
Script to copy just output plots and metadata to a new directory.

TODO:
- option to remove src files after copy?
- handle file globs
'''

import argparse
import os
import shutil


METADATA_DIR = 'metadata'
PLOT_DIR = 'plotOut'

DEFAULT_EXTENSION = '.png'  # For plots
OTHER_EXTENSIONS = ['.tsv']  # For other files


def parse_args():
	'''
	Parses arguments from the command line.

	Returns:
		ArgumentParser namespace: values of variables parsed from the command line
	'''

	parser = argparse.ArgumentParser(description='Recursively copy plots and metadata from one directory to another')

	# Positional arguments
	parser.add_argument('src',
		help='Base directory to copy from')
	parser.add_argument('dest',
		help='Base directory to copy to')

	# Optional arguments
	parser.add_argument('-e', '--extension',
		default=DEFAULT_EXTENSION,
		help='Plot file extension to copy (default: {})'.format(DEFAULT_EXTENSION))
	parser.add_argument('-f', '--force',
		action='store_true',
		help='If set, will copy even if dest directory exists')
	parser.add_argument('-v', '--verbose',
		action='store_true',
		help='If set, prints each copy made')
	parser.add_argument('-d', '--dry-run',
		action='store_true',
		help='If set, performs a dry run printing each copy but does not copy any files')

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_args()

	# Add check to not replace existing files
	if os.path.exists(args.dest) and not args.force and not args.dry_run:
		raise Exception('Destination already exists. Use -f to force a copy.')

	# Prevent replace bug if src not given with directory separator at end
	if not args.src.endswith(os.sep):
		args.src += os.sep

	for root, dirs, files in os.walk(args.src):
		current_dir = os.path.basename(root)
		dest_dirs = root.replace(args.src, '')
		for f in files:
			# Paths for metadata files
			if current_dir == METADATA_DIR:
				pass
			# Paths for plot files
			elif args.extension in f:
				# Bring up nested plot directories (eg plotOut/low_res_plots -> plotOut/)
				while current_dir != PLOT_DIR:
					dest_dirs = os.path.dirname(dest_dirs)
					current_dir = os.path.basename(dest_dirs)

					# If plot directory is not on path then fall back to current directory structure
					if not dest_dirs:
						dest_dirs = root.replace(args.src, '')
						break
			else:
				# Paths for other files
				for e in OTHER_EXTENSIONS:
					if e in f:
						break
				else:
					continue

			src = os.path.join(root, f)
			dest = os.path.join(args.dest, dest_dirs, f)
			status = '{{}} file {} to {}'.format(src, dest)

			# Create directories and copy files
			if args.dry_run:
				print(status.format('Will copy'))
			else:
				dest_dir = os.path.dirname(dest)
				if not os.path.exists(dest_dir):
					os.makedirs(dest_dir)
				shutil.copy2(src, dest)

				if args.verbose:
					print(status.format('Copied'))
