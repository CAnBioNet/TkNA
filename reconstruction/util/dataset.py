from pathlib import Path
import pickle
import shutil
import tempfile
import xarray

class ItemExistsError(Exception):
	pass

class AlreadyLoadedError(Exception):
	pass

class Dataset:
	OBJECTS_FILE_NAME = "objects"

	def __init__(self):
		self.tables = {}
		self.objects = {}
		self.loaded = False

	def get_table_names(self):
		return self.tables.keys()

	def add_table(self, name, table):
		if name in self.tables:
			raise ItemExistsError(f"Table with name \"{name}\" already exists in dataset")
		self.tables[name] = table

	def get_table(self, name, absent_ok=False):
		if name not in self.tables:
			if absent_ok:
				return None
			else:
				raise KeyError(f"Table with name \"{name}\" does not exist in dataset")
		return self.tables[name]

	def add_object(self, key, obj):
		if key in self.objects:
			raise ItemExistsError(f"Object with key \"{key}\" already exists in dataset")
		self.objects[key] = obj

	def get_object(self, key, absent_ok=False):
		if key not in self.objects:
			if absent_ok:
				return None
			else:
				raise KeyError(f"Object with key \"{key}\" does not exist in dataset")
		return self.objects[key]

	def load_from_file(self, file_path):
		if self.loaded:
			raise AlreadyLoadedError("Cannot load to the same dataset more than once")
		temp_dir = tempfile.TemporaryDirectory()
		temp_dir_path = Path(temp_dir.name)
		shutil.unpack_archive(file_path, temp_dir_path, "zip")
		for table_file_path in temp_dir_path.glob("*.cdf"):
			name = table_file_path.stem
			self.tables[name] = xarray.open_dataarray(table_file_path)
		with open(temp_dir_path / self.OBJECTS_FILE_NAME, "rb") as objects_file:
			self.objects = pickle.load(objects_file)
		temp_dir.cleanup()
		self.loaded = True

	def write_to_file(self, file_path, make_parent=False):
		if make_parent:
			file_path.parent.mkdir(parents=True, exist_ok=True)
		temp_dir = tempfile.TemporaryDirectory()
		temp_dir_path = Path(temp_dir.name)
		for name, table in self.tables.items():
			table.to_netcdf(temp_dir_path / f"{name}.cdf")
		with open(temp_dir_path / self.OBJECTS_FILE_NAME, "wb") as objects_file:
			pickle.dump(self.objects, objects_file)
		shutil.make_archive(file_path.parent / file_path.stem, "zip", temp_dir_path)
		temp_dir.cleanup()

