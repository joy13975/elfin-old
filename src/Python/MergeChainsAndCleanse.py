#!/usr/bin/env python

import sys
import argparse
from utils import *
import Bio.PDB

def main():
	ap = argparse.ArgumentParser(description='Merge all chains in a PDB and remove 1H, 2H, 3H, and OXT')

	ap.add_argument('--input')
	ap.add_argument('--output', default='')
	
	if len(sys.argv) == 1:
		ap.print_help()
		sys.exit(1)

	args = ap.parse_args()

	inputExt = args.input[args.input.rfind('.'):]
	die(inputExt != '.pdb', 'Input is not a PDB file')

	# Start the process
	pdb = readPdb('dirtyPdb', args.input)
	model = pdb.child_list[0]
	newChain = Bio.PDB.Chain.Chain('A')

	oldChains = model.child_list

	# Add all residues to a new chain with correct ascending IDs
	rid = 1
	for oc in oldChains:
		for r in oc:
			# Cleanse unwanted atoms
			unwantedAtoms = []
			for a in r.child_list:
				 if a.name == '1H' or a.name == '2H' or a.name == '3H' or a.name == 'OXT':
					unwantedAtoms.append(a.name)

			for ua in unwantedAtoms:
				r.detach_child(ua)

			# Assign new residue ID
			r.id = (r.id[0], rid, r.id[2])
			newChain.add(r)
			rid += 1

	# Remove old chains from model
	ocIds = []
	for oc in oldChains:
		ocIds.append(oc.id)

	for ocId in ocIds:
		model.detach_child(ocId)

	# Add new chain to model
	model.add(newChain)

	# Write to output
	if args.output == '':
		args.output = args.input.replace('.pdb', '_mc.pdb')

	savePdb(pdb, args.output)
	print 'Done and saved to {}'.format(args.output)

if __name__ == '__main__':
	main()