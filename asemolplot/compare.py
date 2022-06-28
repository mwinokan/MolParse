
def compareSystems(system1,system2):
  import mout
  import mcol

  use_cs4qb = True

  try:
    import cs4qb
  except:
    use_cs4qb = False
    mout.warningOut("Not using "+mcol.func+"cs4qb")

  if use_cs4qb:

    import os

    # Get CHARMM36 location
    try:
      c36_path = os.environ["C36CUSTOM"]
    except KeyError:
      mout.errorOut("$C36CUSTOM variable not set",fatal=True)

    # Get FF filepaths
    ffb_filepath = c36_path+"/ffbonded.itp"
    ffnb_filepath = c36_path+"/ffnonbonded.itp"
    rtp_filepath = c36_path+"/merged.rtp"
    hdb_filepath = c36_path+"/merged.hdb"

    # Check if accessible
    for file in (ffb_filepath,rtp_filepath,hdb_filepath):
      if not os.access(file, os.R_OK):
        mout.errorOut(file+" cannot be accessed",fatal=True)

    from cs4qb.gmx_ff import getFFBondList

    system_data = []
    system_data.append({})
    system_data.append({})

    for index,system in enumerate([system1,system2]):

      value_list = []
      bond_list = []
      improper_list = []

      for residue in system.residues:

        this_bond_list, this_improper_list = getFFBondList(residue,rtp_filepath,use_atomtypes=False)

        for bond in this_bond_list:

          if len(bond) == 2:
            bond_list.append([residue.name+"_"+b for b in bond])
            value_list.append(residue.get_distance(bond[0],bond[1]))
          elif len(bond) == 3:
            bond_list.append([residue.name+"_"+b for b in bond])
            value_list.append(residue.get_angle(bond[0],bond[1],bond[2]))

          # else:
            # value_list.append(None)

      system_data[index]['bond_list'] = bond_list
      system_data[index]['values'] = value_list

    import json

    # Write the CJSON files
    with open('system1.json', 'w') as f: json.dump(system_data[0], f,indent=4)
    with open('system2.json', 'w') as f: json.dump(system_data[1], f,indent=4)

    system_data = compare_bond_stat_lists(system_data,dataFile=True)
    
    with open('compare.json', 'w') as f: json.dump(system_data[2], f,indent=4)

def compare_bond_stat_lists(system_data,precision=6,diffPrecision=2,dataFile=False):

  import mout
  import mcol

  diff_list = []
  pcnt_diff_list = []

  if dataFile:
    import os
    datfile = "compare.dat"
    f_dat = open(datfile,'w')
    f_dat.write('BOND, DIFFERENCE, % \n')
    f_dat.close()
    f_dat = open(datfile,'a')
  else:
    f_dat = None

  for index,bond in enumerate(system_data[0]['bond_list']):

    assert system_data[0]['bond_list'][index] == system_data[1]['bond_list'][index]

    if len(bond) == 2:
      bondStr = bond[0]+" - "+bond[1]
    elif len(bond) == 3:
      bondStr = bond[0]+" - "+bond[1]+" - "+bond[2]

    diff, pcnt_diff = mout.differenceOut(bondStr,
                                         system_data[0]['values'][index],
                                         system_data[1]['values'][index],
                                         valCol=mcol.result,
                                         precision=precision,
                                         diffPrecision=diffPrecision,
                                         unit="Angstroms",
                                         dataFile=f_dat)
    
    diff_list.append(diff)
    pcnt_diff_list.append(pcnt_diff)

  diff_avg = 0
  for diff in pcnt_diff_list:
    diff_avg += diff
  diff_avg /= len(pcnt_diff_list)

  if len(system_data) == 2:
    system_data.append({})
  system_data[2]['bond_list'] = system_data[0]['bond_list']
  system_data[2]['diff'] = diff_list
  system_data[2]['pcnt_diff'] = pcnt_diff_list
  system_data[2]['avg_diff'] = diff_avg

  f_dat.write('AVERAGE, N/A, '+str(round(diff_avg,2))+' \n')

  return system_data

  # print(system1.res_names)
  # print(system2.res_names)

  # print(system1.positions)
  # print(system2.positions)

def euclid_dist(system1,system2):

  import numpy as np

  all_positions1=[]
  all_positions2=[]

  for i,atom in enumerate(system1.atoms):
    all_positions1 += atom.position
    all_positions2 += system2.atoms[i].position

  return np.linalg.norm(np.array(all_positions2) - np.array(all_positions1))
