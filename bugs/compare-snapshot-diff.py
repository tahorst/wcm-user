#! /usr/bin/env python

"""
Compares current Jenkins files to snapshot directory on Sherlock to find
differences between successful builds and seg fault builds.
"""

import os
import subprocess


INSTALL_PATH = '/home/groups/mcovert/downloads/jenkins/'
SNAPSHOT_PATH = os.path.join(INSTALL_PATH, '.snapshot', 'groups.monthly.201912')
DIR_SKIP = 'builds'
FILE_SKIP = {
	'hudson.plugins.disk_usage.DiskUsageProjectActionFactory.xml',
	'disk-usage.xml',
	}


if __name__ == '__main__':
	for root, dirs, files in os.walk(INSTALL_PATH):
		# Skip undesired directory if on root path
		if DIR_SKIP in root:
			continue

		# Get same snapshot path as root
		additional_dir = root.split(INSTALL_PATH)[1]
		snapshot = os.path.join(SNAPSHOT_PATH, additional_dir)

		# Check diff for each file found
		for name in files:
			# Skip files that are being modified and cause diff to hang
			if name in FILE_SKIP:
				continue

			file1 = os.path.join(root, name)
			file2 = os.path.join(snapshot, name)

			# Get diff from files
			p = subprocess.Popen(
				['diff', '-yZw', '--suppress-common-lines', file1, file2],
				stdout=subprocess.PIPE
				)
			p.wait()
			diff, _ = p.communicate()
			if diff:
				print('\n*** {} ***'.format(os.path.join(additional_dir, name)))
				print(diff)
