#!/usr/bin/env python3

import mout
import molparse as mp
from mout.testing import TestStatus
from molparse.amino import INTERACTION_SITES, BB_NAMES, RES_NAMES, BB_INTERACTION_SITES

import os
import glob
from pprint import pprint


@mout.debug_log
def check_bb_interactions():
    test_status = TestStatus()

    ### Backbone interactions

    for atom_name in BB_NAMES:

        # get all interactions for those names
        interactions = get_interactions_by_atom_name(atom_name)
        if interactions:
            mout.var(f'interactions({atom_name})', interactions.keys())

        # warn if BB interactions present
        if interactions:
            test_status.warning('BB interactions present, move to AMINO_INTERACTIONS["BB"]')

    return test_status


@mout.debug_log
def check_interaction_source():
    test_status = TestStatus()

    for res_name, interactions in INTERACTION_SITES.items():

        for interaction in interactions:

            if 'source' not in interaction:
                test_status.error(f'No source field in {res_name}: {interaction["type"]} {interaction["atoms"]}')
                continue

            if 'unsure' in interaction['source']:
                test_status.warning(f'unsure source for {res_name}: {interaction["type"]} {interaction["atoms"]}')

    return test_status


@mout.debug_log
def check_interaction_pair(int1, int2):
    test_status = TestStatus()

    # for each interaction in an interaction pair check that there is another interaction with the same atoms

    for res_name, interactions in INTERACTION_SITES.items():

        for interaction in interactions:

            if interaction['type'] == int1:
                matches = [i for i in interactions if i['type'] == int2 and i['atoms'] == interaction['atoms']]
                if len(matches) != 1:
                    test_status.warning(f'{res_name}: No {int2} for {interaction["atoms"]}')

            elif interaction['type'] == int2:
                matches = [i for i in interactions if i['type'] == int1 and i['atoms'] == interaction['atoms']]
                if len(matches) != 1:
                    test_status.warning(f'{res_name}: No {int1} for {interaction["atoms"]}')

    return test_status


@mout.debug_log
def check_pi_stack_pair():
    test_status = TestStatus()

    # for each interaction in an interaction pair check that there is another interaction with the same atoms

    for res_name, interactions in INTERACTION_SITES.items():

        for interaction in interactions:

            # every pi_stack should be able to do a pi_cation
            if interaction['type'] == 'pi_stacking':
                matches = [i for i in interactions if i['type'] == 'pi_cation' and i['atoms'] == interaction['atoms']]
                if len(matches) != 1:
                    test_status.warning(f'{res_name}: No pi_cation for {interaction["atoms"]}')

            # only pi_cation interactions with the ring on the protein infer pi_stacking
            elif interaction['type'] == 'pi_cation':
                if len(interaction['atoms']) < 4:
                    continue
                matches = [i for i in interactions if i['type'] == 'pi_stacking' and i['atoms'] == interaction['atoms']]
                if len(matches) != 1:
                    test_status.warning(f'{res_name}: No pi_stacking for {interaction["atoms"]}')

    return test_status


def get_interactions_by_atom_name(atom_name):
    result = {}

    for res_name, interactions in INTERACTION_SITES.items():

        for interaction in interactions:

            if atom_name in interaction['atoms']:

                int_type = interaction['type']

                if int_type not in result:
                    result[int_type] = []

                if res_name not in result[int_type]:
                    result[int_type].append(res_name)

    return result


@mout.debug_log
def check_amino_interaction_sites_method():
    test_status = TestStatus()

    ref_dir = f'{os.path.dirname(__file__)}/../molparse/ref'

    for res_name in RES_NAMES:

        files = glob.glob(f'{ref_dir}/{res_name}.pdb')

        if len(files) == 0:
            test_status.error(f'no reference file for {res_name}')
            continue

        elif len(files) == 0:
            test_status.error(f'multiple reference files for {res_name}')
            continue

        system = mp.parse(files[0], verbosity=0)

        interactions = system['r0'].interaction_sites

        if not interactions:
            test_status.warning(f'No interactions for {res_name}')

    return test_status


def main():
    mout.header(__file__)

    # used to generate INTERACTION_SITES['BB']
    check_bb_interactions()

    check_interaction_source()
    check_interaction_pair('hydrogen_donor', 'water_donor')
    check_interaction_pair('hydrogen_acceptor', 'water_acceptor')
    check_pi_stack_pair()

    # check that histidine interactions are identical?

    check_amino_interaction_sites_method()

    TestStatus.summary()


if __name__ == '__main__':
    main()
