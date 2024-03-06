#!/usr/bin/env python3

from mout import debug_log
from mout.testing import TestStatus

from mlog import setup_logger
logger = setup_logger('res_type')

from molparse.residue import add_res_type, res_type

@debug_log
def check_res_type(res_name):
	test_status = TestStatus()

	logger.info(f"{res_type(res_name)=}")

	return test_status

def main():

	check_res_type('DG')
	check_res_type('TYR')
	check_res_type('CRO')
	check_res_type('CROW')
	check_res_type('COW')
	check_res_type('DOG')

	add_res_type('COW','PRO')
	add_res_type('CROW','LIP')
	add_res_type('DOG','LIG')
	
	check_res_type('CROW')
	check_res_type('COW')
	check_res_type('DOG')
	
if __name__ == '__main__':
	main()
