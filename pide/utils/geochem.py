import numpy as np
import matplotlib.pyplot as plt
import json

class TASClassifier:
	"""
	Total Alkali-Silica (TAS) diagram classifier for volcanic rocks.
	
	This implementation uses the LeMaitre et al. classification scheme
	
	Field codes and their corresponding rock names (volcanic/plutonic):
	- F: Foidite/Foidolite
	- Pc: Picrite/Peridotgabbro
	- U1: Tephrite/Foid Gabbro
	- U2: Phonotephrite/Foid Monzodiorite
	- U3: Tephriphonolite/Foid Monzosyenite
	- Ba: Alkalic Basalt/Alkalic Gabbro
	- Bs: Subalkalic Basalt/Subalkalic Gabbro
	- O1: Basaltic Andesite/Gabbroic Diorite
	- O2: Andesite/Diorite
	- O3: Dacite/Granodiorite
	- S1: Trachybasalt/Monzogabbro
	- S2: Basaltic Trachyandesite/Monzodiorite
	- S3: Trachyandesite/Monzonite
	- T1: Trachyte/Syenite
	- T2: Trachydacite/Quartz Monzonite
	- R: Rhyolite/Granite
	"""
	
	def __init__(self, config_file=None):
		"""
		Initialize the TAS classifier.
		
		Parameters:
		-----------
		config_file : str, optional
			Path to JSON configuration file. If not provided, uses built-in configuration.
		"""
		if config_file:
			self.load_config(config_file)
		else:
			# Use the built-in configuration based on the provided JSON
			self._load_default_config()
	
	def _load_default_config(self):
		"""Load the default TAS configuration."""
		self.config = {
			"name": "TAS",
			"axes": {
				"x": "SiO2",
				"y": "Na2O + K2O"
			},
			"fields": {
				"Ba": {
					"name": ["Alkalic Basalt", "Alkalic Gabbro"],
					"poly": [[45, 5], [52, 5], [45, 2]]
				},
				"Bs": {
					"name": ["Subalkalic Basalt", "Subalkalic Gabbro"],
					"poly": [[45, 2], [52, 5], [52, 0], [45, 0]]
				},
				"F": {
					"name": ["Foidite", "Foidolite"],
					"poly": [[35, 9], [37, 14], [52.5, 18], [52.5, 14], [48.4, 11.5], 
							[45, 9.4], [41, 7], [41, 3], [37, 3]]
				},
				"O1": {
					"name": ["Basaltic Andesite", "Gabbroic Diorite"],
					"poly": [[52, 0], [52, 5], [57, 5.9], [57, 0]]
				},
				"O2": {
					"name": ["Andesite", "Diorite"],
					"poly": [[57, 0], [57, 5.9], [63, 7], [63, 0]]
				},
				"O3": {
					"name": ["Dacite", "Granodiorite"],
					"poly": [[63, 0], [63, 7], [69, 8], [77.3, 0]]
				},
				"Pc": {
					"name": ["Picrite", "Peridotgabbro"],
					"poly": [[41, 3], [45, 3], [45, 2], [45, 0], [41, 0]]
				},
				"Ph": {
					"name": ["Phonolite", "Foid Syenite"],
					"poly": [[52.5, 14], [52.5, 18], [57, 18], [63, 16.2], 
							[61, 13.5], [57.6, 11.7]]
				},
				"R": {
					"name": ["Rhyolite", "Granite"],
					"poly": [[69, 8], [69, 13], [85.9, 6.8], [87.5, 4.7], [77.3, 0]]
				},
				"S1": {
					"name": ["Trachybasalt", "Monzogabbro"],
					"poly": [[45, 5], [49.4, 7.3], [52, 5]]
				},
				"S2": {
					"name": ["Basaltic Trachyandesite", "Monzodiorite"],
					"poly": [[49.4, 7.3], [53, 9.3], [57, 5.9], [52, 5]]
				},
				"S3": {
					"name": ["Trachyandesite", "Monzonite"],
					"poly": [[53, 9.3], [57.6, 11.7], [61, 8.6], [63, 7], [57, 5.9]]
				},
				"T1": {
					"name": ["Trachyte", "Syenite"],
					"poly": [[57.6, 11.7], [61, 13.5], [63, 16.2], [69, 13], [61, 8.6]]
				},
				"T2": {
					"name": ["Trachydacite", "Quartz Monzonite"],
					"poly": [[61, 8.6], [69, 13], [69, 8], [63, 7]]
				},
				"U1": {
					"name": ["Tephrite", "Foid Gabbro"],
					"poly": [[41, 3], [41, 7], [45, 9.4], [49.4, 7.3], [45, 5], [45, 3]]
				},
				"U2": {
					"name": ["Phonotephrite", "Foid Monzodiorite"],
					"poly": [[45, 9.4], [48.4, 11.5], [53, 9.3], [49.4, 7.3]]
				},
				"U3": {
					"name": ["Tephriphonolite", "Foid Monzosyenite"],
					"poly": [[48.4, 11.5], [52.5, 14], [57.6, 11.7], [53, 9.3]]
				}
			}
		}
	
	def load_config(self, config_file):
		"""
		Load TAS configuration from a JSON file.
		
		Parameters:
		-----------
		config_file : str
			Path to the JSON configuration file
		"""
		with open(config_file, 'r') as f:
			self.config = json.load(f)
	
	def point_in_polygon(self, x, y, polygon):
		"""Check if a point (x, y) is inside a polygon."""
		n = len(polygon)
		inside = False
		
		p1x, p1y = polygon[0]
		for i in range(1, n + 1):
			p2x, p2y = polygon[i % n]
			if y > min(p1y, p2y):
				if y <= max(p1y, p2y):
					if x <= max(p1x, p2x):
						if p1y != p2y:
							xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
						if p1x == p2x or x <= xinters:
							inside = not inside
			p1x, p1y = p2x, p2y
		
		return inside
	
	def classify(self, sio2, na2o, k2o, rock_type='volcanic'):
		"""
		Classify a rock based on its SiO2 and alkali content.
		
		Parameters:
		-----------
		sio2 : float
			SiO2 content in weight %
		na2o : float
			Na2O content in weight %
		k2o : float
			K2O content in weight %
		rock_type : str
			'volcanic' or 'plutonic' to get the appropriate rock name
		
		Returns:
		--------
		str
			Rock classification name
		"""
		total_alkalis = na2o + k2o
		
		# Check each field to see if the point falls within it
		for field_code, field_data in self.config['fields'].items():
			if field_code in ['nan', 'none']:
				continue
			
			polygon = field_data['poly']
			if polygon and self.point_in_polygon(sio2, total_alkalis, polygon):
				names = field_data['name']
				if isinstance(names, list) and len(names) >= 2:
					return names[0] if rock_type == 'volcanic' else names[1]
				elif isinstance(names, list):
					return names[0]
				else:
					return names
		
		# If no match found, check if it's outside the typical range
		if sio2 < 35:
			return "Ultrabasic (SiO2 < 35%)"
		elif sio2 > 87.5:
			return "Ultra-acidic (SiO2 > 87.5%)"
		else:
			return "Unclassified"
	
	def get_field_code(self, sio2, na2o, k2o):
		"""
		Get the field code (abbreviation) for a rock composition.
		
		Parameters:
		-----------
		sio2 : float
			SiO2 content in weight %
		na2o : float
			Na2O content in weight %
		k2o : float
			K2O content in weight %
		
		Returns:
		--------
		str
			Field code (e.g., 'Ba', 'R', 'Ph')
		"""
		total_alkalis = na2o + k2o
		
		# Check each field to see if the point falls within it
		for field_code, field_data in self.config['fields'].items():
			if field_code in ['nan', 'none']:
				continue
				
			polygon = field_data['poly']
			if polygon and self.point_in_polygon(sio2, total_alkalis, polygon):
				return field_code
		
		return "N/A"
	
def classify_tas_diagram(sio2, na2o, k2o):

	"""
	Classify rock samples using the Total Alkali-Silica (TAS) diagram.
	
	Performs TAS classification on multiple rock samples based on their 
	SiO2 and total alkali (Na2O + K2O) content. Returns field codes and 
	both volcanic and plutonic rock names according to the LeMaitre et al. 
	classification scheme.
	
	Parameters
	----------
	sio2 : array_like
		SiO2 content in weight percent for each sample. Must be in the 
		range [30, 90] for meaningful classification.
	na2o : array_like
		Na2O content in weight percent for each sample. Must be non-negative.
	k2o : array_like
		K2O content in weight percent for each sample. Must be non-negative.
		
	Returns
	-------
	dict
		Dictionary containing three keys:
		
		- 'field_codes' : list of str
			TAS field codes for each sample (e.g., 'Ba', 'R', 'Ph').
			Returns 'N/A' if sample falls outside defined fields.
		- 'volcanic_names' : list of str
			Volcanic rock names for each sample (e.g., 'Basalt', 'Rhyolite').
		- 'plutonic_names' : list of str
			Plutonic rock names for each sample (e.g., 'Gabbro', 'Granite').
	
	Examples
	--------
	> sio2 = [48.5, 55.0, 70.0]
	> na2o = [2.3, 3.5, 4.0]
	> k2o = [1.2, 2.0, 4.5]
	> results = classify_tas_diagram(sio2, na2o, k2o)
	> print(results['field_codes'])
	"""

	classifier = TASClassifier()

	field_codes = []
	volcanic_names = []
	plutonic_names = []
	for i in range(len(sio2)):
		fc = classifier.get_field_code(sio2[i], na2o[i], k2o[i])
		vn = classifier.classify(sio2[i], na2o[i], k2o[i], 'volcanic')
		pn = classifier.classify(sio2[i], na2o[i], k2o[i], 'plutonic')
		field_codes.append(fc)
		volcanic_names.append(vn)
		plutonic_names.append(pn)

	return {'field_codes':field_codes, 'volcanic_names':volcanic_names, 'plutonic_names': plutonic_names}