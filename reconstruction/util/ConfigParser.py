import json

class ExpectedValueError(Exception):
	def __init__(self, key):
		message = "Expected a value for \"{}\", but no value was given.".format(key)
		super().__init__(message)

class InvalidTypeError(Exception):
	def __init__(self, key, expectedType, recievedType, typeList=False):
		if typeList:
			message = "Invalid type for \"{}\". Expected one of \"{}\", but got \"{}\".".format(key, expectedType, recievedType)
		else:
			message = "Invalid type for \"{}\". Expected \"{}\", but got \"{}\".".format(key, expectedType, recievedType)
		super().__init__(message)

class ValidationFailedError(Exception):
	pass

class NoDefault:
	pass

class UseSubDefaults:
	pass

class ConfigItem:
	def __init__(self, key, default=NoDefault()):
		self.key = key
		self.default = default

	def getDefault(self):
		if isinstance(self.default, NoDefault):
			raise ExpectedValueError(self.key)
		return self.default

	def _parse(self, jsonValue):
		raise NotImplementedError()

	def parse(self, jsonConfig):
		if self.key in jsonConfig:
			return self._parse(jsonConfig[self.key])
		return self.getDefault()

class MultiTypeConfigItem(ConfigItem):
	def __init__(self, key, possibleTypes, default=NoDefault()):
		super().__init__(key, default)
		self.possibleTypes = possibleTypes

	def _parse(self, jsonValue):
		for possibleType in self.possibleTypes:
			if isinstance(jsonValue, possibleType):
				return jsonValue
		raise InvalidTypeError(self.key, self.possibleTypes, type(jsonValue), typeList=True)

class TypedConfigItem(MultiTypeConfigItem):
	def __init__(self, key, expectedType, default=NoDefault()):
		super().__init__(key, [expectedType], default)

class TypeOrSubconfigConfigItem(ConfigItem):
	def __init__(self, key, type_, subconfig, default=NoDefault()):
		super().__init__(key, default)
		self.type_ = type_
		self.subconfig = subconfig

	def getDefault(self):
		if isinstance(self.default, UseSubDefaults):
			return self.subconfig.getDefault()
		return super().getDefault()

	def _parse(self, jsonValue):
		if isinstance(jsonValue, self.type_):
			return jsonValue
		if isinstance(jsonValue, dict):
			return self.subconfig.parse(jsonValue)
		raise InvalidTypeError(self.key, self.type_, type(jsonValue))

class Config:
	def __init__(self, configItems, validator=None):
		self.configItems = configItems
		self.validator = validator

	def getDefault(self):
		defaults = {}
		for item in self.configItems:
			 defaults[item.key] = item.getDefault()
		return defaults

	def parse(self, jsonConfig):
		values = {}
		for item in self.configItems:
			values[item.key] = item.parse(jsonConfig)
		return values

	def validate(self, config):
		if self.validator is not None:
			return self.validator(config)
		else:
			return config

def parseConfig(configSpec, jsonConfig):
	config = configSpec.parse(jsonConfig)
	config = configSpec.validate(config)
	return config

def parseConfigFile(configSpec, configFilePath):
	configFile = open(configFilePath, "r")
	jsonConfig = json.load(configFile)
	configFile.close()

	return parseConfig(configSpec, jsonConfig)
	
