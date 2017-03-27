#!/usr/bin/env python

from __future__ import print_function
import os

import argparse
from argparse import RawTextHelpFormatter


def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None

def main():
	parser = argparse.ArgumentParser(
			formatter_class=RawTextHelpFormatter,
			description='Check dependencies')
	parser.add_argument('deps', metavar='PROG', nargs='*', help='programs to check')

	args = parser.parse_args()

	print('Checking dependencies:')
	for dep in args.deps:
		if which(dep):
			print(' ........ '.join([dep, 'Found {}'.format(which(dep)), '[OK]']))
		else:
			print(' ........ '.join([dep, 'Found {}'.format(which(dep)), '[ERROR]']))

if __name__ == "__main__":
	main()