#!/usr/bin/env python

import argparse, sys
# from utils import *

def main():
	ap = argparse.ArgumentParser(description='Template Python script');
	# ap.add_argument('--input')

	# Exit if no argument given
	if len(sys.argv) == 1:
		ap.print_help()
		sys.exit(1)

	args = ap.parse_args()

	# Do stuff...

if __name__ == '__main__':
	main()