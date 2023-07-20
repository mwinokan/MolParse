
from mwin.curses import CursesApp, Button, Text

import curses

from .group import AtomGroup
from .system import System
from .chain import Chain
from .residue import Residue
from .atom import Atom
from .amino import AminoAcid
from .nucleic import NucleicAcid

""" To Do's

* Simplify chains of single atom/residue types repeated over and over
* Distinguish between hetatm/backbone/sidechain atoms
* Add more capabilities to context windows
* Add cursor selection possibility

"""

def tree(obj):
	"""Open a CLI App showing the hierarchical tree of a MolParse System."""
	
	viewer = TreeViewer(obj)

class TreeViewer(CursesApp):
	"""CLI App showing the hierarchical tree of a MolParse System. Use molparse.tree(sys) to start."""
	
	def __init__(self, obj):
		
		super(TreeViewer, self).__init__(debug=False)

		if isinstance(obj,str):
			from .io import parse
			obj = parse(obj)

		self._obj = obj

		if isinstance(obj,System):

			obj._expand = True

			if obj.num_chains == 1:
				obj.chains[0]._expand = True

				if obj.chains[0].num_residues == 1:
					obj.chains[0].residues[0]._expand = True

		if isinstance(obj, AtomGroup):
			obj._expand = True

		else:
			import mout
			mout.errorOut(f"TreeViewer not supported with {self._obj_type}")

		try:
			
			self.drawtree()
			self.firstdraw()

			while True:

				redraw = self.draw()

				if redraw:
					self.clear_drawables()
					self.drawtree()
					self.drawcore()
				
		except KeyboardInterrupt:
			pass
		except curses.error as e:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to a curses error.")
			if 'prefresh() returned ERR' in str(e):
				mout.error('Maybe you resized the terminal window?')
			mout.error(e)
		except KeyError as e:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to a KeyError (maybe plot3d?).")
			mout.error(e)
		except AttributeError as e:
			self.close()
			print()
			import mout
			mout.error("TreeViewer exited due to an AttributeError (maybe outdated MPyTools?).")
			mout.error(e)
		except Exception as e:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to an error.")
			mout.error(e)

		self.close()

	@property
	def obj(self):
		return self._obj

	def drawtree(self):
		max_index = self.obj.children[-1].index
		if self.type_str(self.obj.children[-1]) not in ['System','AtomGroup','Chain']:
			max_number = self.obj.children[-1].number
		else:
			max_number = None
		self.recursive_tree(self.obj,0,max_index,max_number,0)
			
	def recursive_tree(self,parent,line,max_index,max_number,depth=0):

		line = self.object_line(parent, line, max_index, max_number, depth)
		
		if parent.children and parent._expand:
			max_index = parent.children[-1].index
			for child in parent.children:
				line = self.recursive_tree(child,line,max_index,max_number,depth+1,)

		return line

	def has_pad_space(self,num_children):
		if self.pad_h - self.max_padline < num_children:
			return False
		return True

	def spawn_context(self,obj,line,col):

		col = self._last_pressed.col

		items,clickables = obj._context_info

		widgets = []

		for i,(key,prop) in enumerate(items.items()):

			label = Text(f"{key}: ",line+1+i,col,color_pair=self.WHITE_INV|curses.A_BOLD)
			data = Text(f"{prop} ",line+1+i,label.endcol,color_pair=self.WHITE_INV)
			widgets.append(label)
			widgets.append(data)

		for j,(key,func) in enumerate(clickables.items()):

			label = Button(
				app=self,
				line=line+1+i+j,
				col=col,
				active=False,
				enabler=func,
				target=obj,
				disabler=func,
				color_inactive=self.GREEN_INV,
				color_active=self.GREEN_INV,
				name=key)
			widgets.append(label)

		self.context_menu(line,col,widgets)

	def object_line(self,obj,line,max_index,max_number,depth):

		col = 4*depth

		obj._tree = self
		obj._gui_line = line

		# add expand/collapse button
		if not isinstance(obj,Atom) and self.has_pad_space(len(obj.children)):

			f1 = lambda x: x.expand()
			f2 = lambda x: x.collapse()

			button = Button(app=self,
							line=line,
							col=col,
							active=obj._expand,
							enabler=f1,
							color_active=self.RED_INV,
							disabler=f2,
							target=obj,
							name=f' > ',
							activename=f' v ')
			self.add_button(button)
			col = button.endcol + 1

		# no expand functionality
		else:

			text = Text(f'   ',line,col,bold=False)
			self.add_text(text)
			col = text.endcol + 1

		# object type
		text = Text(f'{self.type_str(obj)} ',line,col,bold=True)
		self.add_text(text)

		# object index
		if self.type_str(obj) not in ['System','AtomGroup']:

			# index
			width = len(str(max_index))
			text = Text(f'{str(obj.index).rjust(width)} ',line,text.endcol,color_pair=self.GREEN,bold=True)
			self.add_text(text)

		if self.type_str(obj) not in ['System','AtomGroup','Chain']:
			# number
			width = len(str(max_number))
			text = Text(f'{str(obj.number).rjust(width)} ',line,text.endcol,color_pair=self.GREEN,bold=True)
			self.add_text(text)

		# object name is also a button:
		if self.type_str(obj) not in ['Atom']:
			name = obj.name

			if len(name) > 23:
				name = f'{name[0:20]}...'

		else:
			name = obj.symbol

		f1 = lambda x: x.spawn_context(line,int(text.endcol)-len(name))
		f2 = lambda x: x.hide_context()

		if not isinstance(obj,Atom):
			# obj._show_context = False

			text = Button(app=self,
							line=line,
							col=text.endcol,
							active=obj._show_context,
							enabler=f1,
							color_active=self.CYAN_INV,
							color_inactive=self.CYAN_INV,
							# color_active=self.WHITE_INV|curses.A_BOLD,
							# color_inactive=curses.A_UNDERLINE|self.CYAN|curses.A_BOLD,
							disabler=f2,
							target=obj,
							name=name)
							# activename=f'!{name}!')
			self.add_button(text)

		else:
			text = Text(name,line,text.endcol,color_pair=self.CYAN,bold=True)
			self.add_text(text)

		# residue/chain type
		if self.type_str(obj) not in ['System','Atom','AtomGroup']:
			text = Text(f'(',line,text.endcol+1,bold=True)
			self.add_text(text)
			text = Text(f'{obj.type}',line,text.endcol,color_pair=self.MAGENTA,bold=True)
			self.add_text(text)
			text = Text(f')',line,text.endcol,bold=True)
			self.add_text(text)

		# position
		if self.type_str(obj) in ['Atom']:

			text = Text(f' [{obj.x:7.3f}, {obj.y:7.3f}, {obj.z:7.3f}]',line,text.endcol)
			self.add_text(text)

			text = Text(f' {obj.name}',line,text.endcol,color_pair=self.MAGENTA,bold=True)
			self.add_text(text)

		# extra info
		if self.scr_w >= 80:
			col = self.extra_info(obj, line, text.endcol+1)
		else:
			col = text.endcol

		# plotting buttons
		if self.type_str(obj) in ['System', 'Chain', 'Residue', 'AminoAcid', 'NucleicAcid', 'AtomGroup']:

			col += 1

			f1 = lambda x: x.plot3d()

			button = Button(app=self,
							line=line,
							col=col,
							active=False,
							color_inactive=curses.A_UNDERLINE|curses.A_BOLD,
							enabler=f1,
							disabler=f1,
							target=obj,
							name=f'plotly')
			
			self.add_button(button)
			col = button.endcol + 1
			
			f1 = lambda x: x.view()

			button = Button(app=self,
							line=line,
							col=col,
							active=False,
							color_inactive=curses.A_UNDERLINE|curses.A_BOLD,
							enabler=f1,
							disabler=f1,
							target=obj,
							name=f'ASE')

			self.add_button(button)
			col = button.endcol + 1

		if self.type_str(obj) in ['Residue', 'AminoAcid', 'NucleicAcid', 'Chain']:

			file_out = f'{self.type_str(obj)}_{obj.name}_{obj.index}.pdb'

			f1 = lambda x: x.write(file_out)

			button = Button(app=self,
							line=line,
							col=col,
							active=False,
							color_inactive=curses.A_UNDERLINE|curses.A_BOLD,
							enabler=f1,
							disabler=f1,
							target=obj,
							name=f'PDB')

			self.add_button(button)
			col = button.endcol + 1


		return line + 1

	def extra_info(self,obj,line,col):
		if isinstance(obj,System):
			text = Text(f'contains: {self.quantity_str(obj.num_chains,"chain")}, {self.quantity_str(obj.num_residues,"residue")}, {self.quantity_str(obj.num_atoms,"atom")} ',line,col)
			self.add_text(text)
		elif isinstance(obj,Chain):
			text = Text(f'contains: {self.quantity_str(obj.num_residues,"residue")}, {self.quantity_str(obj.num_atoms,"atom")} ',line,col)
			self.add_text(text)
		elif isinstance(obj,Residue):
			text = Text(f'contains: {self.quantity_str(obj.num_atoms,"atom")} ',line,col)
			self.add_text(text)
		else:
			return col
		return text.endcol

	def quantity_str(self,number,name):
		if number > 1:
			return f'{number} {name}s'
		else:
			return f'{number} {name}'

	def type_str(self,obj):
		if isinstance(obj,System): return "System"
		if isinstance(obj,Chain): return "Chain"
		if isinstance(obj,AminoAcid): return "AminoAcid"
		if isinstance(obj,NucleicAcid): return "NucleicAcid"
		if isinstance(obj,Residue): return "Residue"
		if isinstance(obj,Atom): return "Atom"
		if isinstance(obj,AtomGroup): return "AtomGroup"
