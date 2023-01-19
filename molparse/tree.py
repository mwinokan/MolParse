
from mwin.curses import CursesApp, Button, Text

import curses

from .system import System
from .chain import Chain
from .residue import Residue
from .atom import Atom
from .amino import AminoAcid
from .nucleic import NucleicAcid

def tree(obj):
	"""Open a CLI App showing the hierarchical tree of a MolParse System."""
	
	viewer = TreeViewer(obj)

class TreeViewer(CursesApp):
	"""CLI App showing the hierarchical tree of a MolParse System. Use molparse.tree(sys) to start."""
	
	def __init__(self, obj):
		
		super(TreeViewer, self).__init__()

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

		else:
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
		except curses.error:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to a curses error.")
		except KeyError:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to a KeyError (maybe plot3d?).")
		except AttributeError:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to an AttributeError (maybe outdated MPyTools?).")
		except:
			self.close()
			print()
			import mout
			mout.errorOut("TreeViewer exited due to an error.")

		self.close()

	@property
	def obj(self):
		return self._obj

	def drawtree(self):
		max_index = self.obj.children[-1].index
		self.recursive_tree(self.obj,0,max_index,0)
			
	def recursive_tree(self,parent,line,max_index,depth=0):

		line = self.object_line(parent, line, max_index, depth)
		
		if parent.children and parent._expand:
			max_index = parent.children[-1].index
			for child in parent.children:
				line = self.recursive_tree(child,line,max_index,depth+1,)

		return line

	def object_line(self,obj,line,max_index,depth):

		col = 4*depth

		if not isinstance(obj,Atom):

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

		text = Text(f'{self.type_str(obj)} ',line,col,bold=True)
		self.add_text(text)

		if not isinstance(obj,System):
			width = len(str(max_index))
			text = Text(f'{str(obj.index).rjust(width)} ',line,text.endcol,color_pair=self.GREEN,bold=True)
			self.add_text(text)

		if self.type_str(obj) not in ['Atom']:
			text = Text(f'{obj.name} ',line,text.endcol,color_pair=self.CYAN,bold=True)
			self.add_text(text)
		else:
			text = Text(f'{obj.symbol} ',line,text.endcol,color_pair=self.CYAN,bold=True)
			self.add_text(text)

		if self.type_str(obj) not in ['System','Atom']:
			text = Text(f'(',line,text.endcol,bold=True)
			self.add_text(text)
			text = Text(f'{obj.type}',line,text.endcol,color_pair=self.MAGENTA,bold=True)
			self.add_text(text)
			text = Text(f') ',line,text.endcol,bold=True)
			self.add_text(text)

		if self.type_str(obj) in ['Atom']:

			text = Text(f'[{obj.x:7.3f}, {obj.y:7.3f}, {obj.z:7.3f}]',line,text.endcol)
			self.add_text(text)

			text = Text(f' {obj.name}',line,text.endcol,color_pair=self.MAGENTA,bold=True)
			self.add_text(text)

		if self.scr_w >= 80:
			col = self.extra_info(obj, line, text.endcol)
		else:
			col = text.endcol

		if self.type_str(obj) in ['System','Residue', 'AminoAcid', 'NucleicAcid']:

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
