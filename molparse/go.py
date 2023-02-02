
def plot3d(atoms,extra=[],bonds=[],alpha=1.0):
	"""Render the atoms with plotly graph objects. 
	extra can contain pairs of coordinates to be shown as vectors."""

	import plotly.graph_objects as go
	from ase.data import vdw_radii, atomic_numbers
	from ase.data.colors import jmol_colors

	species = [a.symbol for a in atoms]
	species = list(set(species))

	fig = go.Figure()

	if bonds:

		x = []
		y = []
		z = []

		for a,b in bonds:
			x.append(a[0])
			x.append(b[0])
			x.append(None)
			y.append(a[1])
			y.append(b[1])
			y.append(None)
			z.append(a[2])
			z.append(b[2])
			z.append(None)

		trace = go.Scatter3d(x=x,y=y,z=z,mode='lines',name='bonds',line=dict(color='black',width=16))

		fig.add_trace(trace)

	for s in species:

		atom_subset = [a for a in atoms if a.symbol == s]

		data = {}

		positions = [a.np_pos for a in atom_subset]
		x = [p[0] for p in positions]
		y = [p[1] for p in positions]
		z = [p[2] for p in positions]
		
		data['x'] = x
		data['y'] = y
		data['z'] = z

		atomic_number = atomic_numbers[s]
		size = vdw_radii[atomic_number] * 15
		color = jmol_colors[atomic_number]
		color = (color[0]*256,color[1]*256,color[2]*256,alpha)
		
		data['index'] = [a.index for a in atom_subset]
		data['residue'] = [a.residue for a in atom_subset]
		data['res_number'] = [a.res_number for a in atom_subset]
		data['name'] = [a.name for a in atom_subset]

		customdata = []

		for a in atom_subset:
			customstr = f'name={a.name}<br>index={a.index}<br>residue={a.residue}<br>res_number={a.res_number}'
			customdata.append(customstr)

		trace = go.Scatter3d(x=x,y=y,z=z,mode='markers',name=s,marker=dict(size=size,color=f'rgba{color}',line=dict(color='black',width=2)),customdata=customdata,hovertemplate="%{customdata}<extra></extra>")

		fig.add_trace(trace)

	for i,(a,b) in enumerate(extra):

		trace = go.Scatter3d(x=[a[0],b[0]],y=[a[1],b[1]],z=[a[2],b[2]],name=f'extra[{i}]')

		fig.add_trace(trace)

	fig.update_layout(scene_aspectmode='data')

	fig.show()
