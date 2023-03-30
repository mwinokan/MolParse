
from collections import UserList

class NamedList(UserList):
	
	"""Subclass of collections.UserList which allows for subscription using object names"""

	def __init__(self, inherited_list=[]):
		super(NamedList, self).__init__(inherited_list)

	def __getitem__(self,key):

		# if key is a slice
		if isinstance(key,slice):
			return NamedList(self.data[key])

		# try the key as an integer index
		try:
			key = int(key)
			return self.data[key]

		except ValueError:

			# key may contain space-separated values
			split_key = key.split(" ")

			# if just a word return all matches
			if len(split_key) == 1:
				return [obj for obj in self.data if obj.name == key]

			# if two words
			elif len(split_key) == 2:

				# first word should be a name
				name = split_key[0]

				# second word should be an index
				try:
					index = int(split_key[1])
				except ValueError:
					import mout
					mout.errorOut(f"Could not convert {split_key[1]} to int")
					mout.errorOut("Lookup key should be in the form 'NAME INDEX' for a specific lookup")
					return None

				# find matches
				matches = [obj for i,obj in enumerate(self.data) if obj.name == name and i == index]
				
				if not matches:
					return None
				elif len(matches) > 1:
					mout.errorOut("Multiple matches to a supposedly unique object!",fatal=True)
				else:
					return matches[0]

			else:
				import mout
				mout.errorOut("key must be an integer, object name or 'obj_name index'")
				return None
