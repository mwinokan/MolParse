#!/usr/bin/env python3

import mout
from molparse.io import parse

from mlog import setup_logger
logger = setup_logger('del_res')

from mout.testing import TestStatus

@mout.debug_log
def check_res_delete_names(system, names):
	test_status = TestStatus()

	system.remove_residues(names=names)

	sys_names = system.res_names
	for name in names:
		if name in sys_names:
			test_status.error(f'Residue(s) with {name=} not deleted')

	return test_status

@mout.debug_log
def check_res_delete_indices(system, indices):
	test_status = TestStatus()

	system.remove_residues(indices=indices)

	sys_indices = system.res_indices
	for index in indices:
		if index in sys_indices:
			test_status.error(f'Residue(s) with {index=} not deleted')

	return test_status

@mout.debug_log
def check_res_delete_numbers(system, numbers):
	test_status = TestStatus()

	system.remove_residues(numbers=numbers)

	sys_numbers = system.res_numbers
	for number in numbers:
		if number in sys_numbers:
			test_status.error(f'Residue(s) with {number=} not deleted')

	return test_status

@mout.debug_log
def check_atom_delete_names(system, names):
	test_status = TestStatus()

	system.remove_atoms(names=names)

	sys_names = system.atom_names
	for name in names:
		if name in sys_names:
			test_status.error(f'Residue(s) with {name=} not deleted')

	return test_status

@mout.debug_log
def check_atom_delete_indices(system, indices):
	test_status = TestStatus()

	system.remove_atoms(indices=indices)

	sys_indices = system.atom_indices
	for index in indices:
		if index in sys_indices:
			test_status.error(f'Residue(s) with {index=} not deleted')

	return test_status

@mout.debug_log
def check_atom_delete_numbers(system, numbers):
	test_status = TestStatus()

	system.remove_atoms(numbers=numbers)

	sys_numbers = system.atom_numbers
	for number in numbers:
		if number in sys_numbers:
			test_status.error(f'Residue(s) with {number=} not deleted')

	return test_status

def main():

	# get a test system
	import urllib.request
	entry = '3PJR'

	logger.info(f'Retrieving {entry} from protein data bank...')

	urllib.request.urlretrieve(f'http://files.rcsb.org/download/{entry}.pdb', f'{entry}.pdb')

	sys = parse(f'{entry}.pdb')

	sys.summary()

	check_res_delete_names(sys, ["DT3", "DG", "TYR"])
	check_res_delete_indices(sys, [1, 3, 50, 300, 1000, 2000])
	check_res_delete_numbers(sys, [1, 3, 50, 300, 1000, 2000])

	check_atom_delete_names(sys, ["H41", "N9", "CG3"])
	check_atom_delete_indices(sys, [1, 3, 50, 300, 1000, 2000, 8000, 20000])
	check_atom_delete_numbers(sys, [1, 3, 50, 300, 1000, 2000, 8000, 20000])
	
if __name__ == '__main__':
	main()
