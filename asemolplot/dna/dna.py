
import mout
import mcol

def get_hbond_pairs(residue1,residue2):
	legal = check_legality(residue1,residue2)

	if not legal:
		mout.errorOut("Incompatible DNA residues "+str(residue1)+" & "+str(residue2),fatal=True)

	if residue1.name == "DA":
		DA_N1 = residue1.get_atom(name="N1")
		DA_N6 = residue1.get_atom(name="N6")
		DA_H61 = residue1.get_atom(name="H61")
		DA_H62 = residue1.get_atom(name="H62")
		DT_O4 = residue2.get_atom(name="O4")
		DT_N3 = residue2.get_atom(name="N3")
		DT_H3 = residue2.get_atom(name="H3")
		return [[ DA_N1, DT_N3, DT_H3 ],
				[ DT_O4, DA_N6, DA_H61, DA_H62 ]]
	elif residue1.name == "DT":
		DA_N1 = residue2.get_atom(name="N1")
		DA_N6 = residue2.get_atom(name="N6")
		DA_H61 = residue2.get_atom(name="H61")
		DA_H62 = residue2.get_atom(name="H62")
		DT_O4 = residue1.get_atom(name="O4")
		DT_N3 = residue1.get_atom(name="N3")
		DT_H3 = residue1.get_atom(name="H3")
		return [[ DA_N1, DT_N3, DT_H3 ],
				[ DT_O4, DA_N6, DA_H61, DA_H62 ]]
	elif residue1.name == "DC":
		DC_O2 = residue1.get_atom(name="O2")
		DC_N3 = residue1.get_atom(name="N3")
		DC_N4 = residue1.get_atom(name="N4")
		DC_H41 = residue1.get_atom(name="H41")
		DC_H42 = residue1.get_atom(name="H42")
		DG_N1 = residue2.get_atom(name="N1")
		DG_H1 = residue2.get_atom(name="H1")
		DG_N2 = residue2.get_atom(name="N2")
		DG_H21 = residue2.get_atom(name="H21")
		DG_H22 = residue2.get_atom(name="H22")
		DG_C6 = residue2.get_atom(name="C6")
		DG_O6 = residue2.get_atom(name="O6")
		return [[ DC_O2, DG_N2, DG_H21, DG_H22 ],
				[ DC_N3, DG_N1, DG_H1 ],
				[ DG_O6, DC_N4, DC_H41, DC_H42 ]]
	elif residue1.name == "DG":
		DC_O2 = residue2.get_atom(name="O2")
		DC_N3 = residue2.get_atom(name="N3")
		DC_N4 = residue2.get_atom(name="N4")
		DC_H41 = residue2.get_atom(name="H41")
		DC_H42 = residue2.get_atom(name="H42")
		DG_N1 = residue1.get_atom(name="N1")
		DG_H1 = residue1.get_atom(name="H1")
		DG_N2 = residue1.get_atom(name="N2")
		DG_H21 = residue1.get_atom(name="H21")
		DG_H22 = residue1.get_atom(name="H22")
		DG_C6 = residue1.get_atom(name="C6")
		DG_O6 = residue1.get_atom(name="O6")
		return [[ DC_O2, DG_N2, DG_H21, DG_H22 ],
				[ DC_N3, DG_N1, DG_H1 ],
				[ DG_O6, DC_N4, DC_H41, DC_H42 ]]

	return False

def check_legality(residue1,residue2):

	if residue1.name == "DA":
		if residue2.name == "DT":
			return True
	elif residue1.name == "DT":
		if residue2.name == "DA":
			return True
	elif residue1.name == "DG":
		if residue2.name == "DC":
			return True
	elif residue1.name == "DC":
		if residue2.name == "DG":
			return True

	return False

def get_hbond_indices(residue):

	if residue.name == "DA":
		N1 = residue.get_atom(name="N1")
		N6 = residue.get_atom(name="N6")
		N6 = residue.get_atom(name="N6")