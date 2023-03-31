
from collections import UserList

class NamedList(UserList):
	
	"""Subclass of collections.UserList which allows for subscription using object names and numbers according to

	key: int 		returns list[int]
	key: 'n'+int 	returns list[obj.number == int]
	key: str 		returns list[obj.name == str]
	key: str' 'int  returns list[obj.name == str and obj.index == int]
	key: str' i'int  returns list[obj.name == str and obj.index == int]
	
	"""

	def __init__(self, inherited_list=[]):
		super(NamedList, self).__init__(inherited_list)

	def get_matches(self,name=None,number=None):

		assert any([name,number is not None])

		# name only
		if name and number is not None:
			matches = [obj for obj in self.data if obj.name == name]
			return matches
		
		# name and number
		elif name and number is not None:
			matches = [obj for obj in self.data if obj.number == number and obj.name == name]

			if not matches:
				return None
			elif len(matches) > 1:
				mout.errorOut("Multiple matches to a supposedly unique object!",fatal=True)
			else:
				return matches[0]

		# number only
		elif not name and number is not None:
			matches = [obj for obj in self.data if obj.number == number]

			if not matches:
				return None
			elif len(matches) > 1:
				mout.errorOut("Multiple matches to a supposedly unique object!",fatal=True)
			else:
				return matches[0]

	def __getitem__(self,key):

		# if key is a slice
		if isinstance(key,slice):
			return NamedList(self.data[key])

		# try the key as an integer index
		try:
			key = int(key)
			return self.data[key]

		except ValueError:

			# use 
			if key.startswith('n'):

				try:
					number = int(key[1:])
				except ValueError:
					import mout
					mout.errorOut(f"Could not convert {key[1:]} to int")
					mout.errorOut("Lookup key should be in the form 'n'INDEX for a specific lookup")
					return None

				# find matches
				matches = self.get_matches(number=number)
				return matches

			# key may contain space-separated values
			split_key = key.split(" ")

			# if just a word return all matches
			if len(split_key) == 1:
				return self.get_matches(name=split_key[0])

			# if two words
			elif len(split_key) == 2:

				# first word should be a name
				name = split_key[0]

				if split_key[1].startswith('n'):
					try:
						number = int(split_key[1][1:])
					except ValueError:
						import mout
						mout.errorOut(f"Could not convert {key[1:]} to int")
						mout.errorOut("Lookup key should be in the form NAME' n'INDEX for a specific lookup")
						return None

					# find matches
					matches = self.get_matches(name=name,number=number)
					return matches

				else:

					# second word should be an index
					try:
						index = int(split_key[1])
					except ValueError:
						import mout
						mout.errorOut(f"Could not convert {split_key[1]} to int")
						mout.errorOut("Lookup key should be in the form NAME' 'INDEX for a specific lookup")
						return None

				match = self.data[index]

				assert match.name == name

				return match

			else:
				import mout
				mout.errorOut("key must be an integer, object name or 'obj_name index'")
				return None