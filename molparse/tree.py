
from mwin.curses import CursesApp, Button, Text

import curses

from .system import System
from .chain import Chain
from .residue import Residue
from .atom import Atom

def tree(obj):
	
	viewer = TreeViewer(obj)

class TreeViewer(CursesApp):
	"""CLI View of a MolParse object"""
	
	def __init__(self, obj):
		
		super(TreeViewer, self).__init__()

		self._obj = obj

		if isinstance(obj,System):

			obj._expand = True

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
			mout.errorOut("TreeViewer exited due to a curses error.")

		self.close()

	@property
	def obj(self):
		return self._obj

	def drawtree(self):
		self.recursive_tree(self.obj,0,0)
			
	def recursive_tree(self,parent,line,depth=0):

		line = self.object_line(parent, line, depth)
		
		if parent.children and parent._expand:
			for child in parent.children:
				line = self.recursive_tree(child,line,depth+1)

		return line

	def object_line(self,obj,line,depth):

		text = Text(f'{self.type_str(obj)}: ',line,depth*2)
		self.add_text(text)

		if not isinstance(obj,Atom):

			f1 = lambda x: x.expand()
			f2 = lambda x: x.collapse()

			button = Button(app=self,
							name=f' {obj} ',
							line=line,
							col=text.endcol,
							active=obj._expand,
							enabler=f1,
							disabler=f2,
							target=obj)
			self.add_button(button)
			col = button.endcol

		else:

			text = Text(name=f'{obj}',
						line=line,
						col=text.endcol)
			self.add_text(text)
			col = text.endcol

		self.extra_info(obj, line, col)

		return line + 1

	def extra_info(self,obj,line,col):
		if isinstance(obj,System):
			text = Text(f' #chains={obj.num_chains}, #residues={obj.num_residues}, #atoms={obj.num_atoms}',line,col)
			self.add_text(text)
		elif isinstance(obj,Chain):
			text = Text(f' #residues={obj.num_residues}, #atoms={obj.num_atoms}',line,col)
			self.add_text(text)
		elif isinstance(obj,Residue):
			text = Text(f' #atoms={obj.num_atoms}',line,col)
			self.add_text(text)

	def type_str(self,obj):
		if isinstance(obj,System): return "System"
		if isinstance(obj,Chain): return "Chain"
		if isinstance(obj,Residue): return "Residue"
		if isinstance(obj,Atom): return "Atom"
