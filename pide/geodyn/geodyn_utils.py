#!/usr/bin/env python3

def rewrite_uw_xmf_file(input_xmf, iter_index, output_dir, field_strs, field_search = 'projDensityField', dtype = 'Cell'):

	"""A function to add new fields to Underworld xmf file format.
	Input:
	str: input_xmf - filename for the input xmf file.
	int: iter_index - iteration number of the xmf file.
	str: output_dir - output filename
	list: field_strs - list of strings where the newly added fields will be named.
	str: field_search - name of the field that the new fields will be added after in the tree.
	str: dtype - type of the field to be added to the xmf file. || 'Cell' or 'Node'.
	
	Output:
	Writes out edited xmf file at output_dir.
	"""
	
	import os
	import xml.etree.ElementTree as ET
	import xml.dom.minidom as minidom

	if type(field_strs) is not list:
		raise KeyError('field_strs parameter has to be entered as a list of strings that would be added to the xmf file as a new parameter. Even if you have one parameter to add, insert it as ["Field_str1"]')

	# Load the XMF file
	tree = ET.parse(input_xmf)
	root = tree.getroot()
	
	# Find the DeformCondField element and its parent
	field_search_element = None
	parent_element = None
	
	for parent in root.iter():
		for child in parent:
			if child.attrib.get('Name') == field_search:
				field_search_element = child
				parent_element = parent
				break
		if field_search_element is not None:
			break
			
	for parent in root.iter():
		for child in parent:
			if child.attrib.get('Center') == dtype:
				dit = child.find('DataItem')
				if dit is not None:
					dimensions = dit.get('Dimensions')
					break
					
	if field_search_element is not None and parent_element is not None:
		# Create the new ABCondField element
		for i in range(0,len(field_strs)):
			new_field = ET.Element('Attribute', attrib={
				'Type': 'Scalar',
				'Center': dtype,
				'Name': field_strs[i]
			})
		
			# Add a DataItem child to ABCondField
			data_item = ET.SubElement(new_field, 'DataItem', attrib={
				'Format': 'HDF',
				'NumberType': 'Float',
				'Precision': '8',
				'Dimensions': dimensions
			})
			data_item.text = field_strs[i] +'-' + str(iter_index) + '.h5:/data'
		
			# Insert ABCondField after DeformCondField
			index = list(parent_element).index(field_search_element)
			parent_element.insert(index + i + 1, new_field)
			
			output_name = os.path.join(os.getcwd(),output_dir,input_xmf)
			
			# Convert the tree to a string
			xml_str = ET.tostring(root, encoding='unicode')
	
			# Pretty-print the XML
			pretty_xml_str = minidom.parseString(xml_str).toprettyxml(indent="    ")
		
			# Save the pretty-printed XML to a file
			with open(output_name, 'w') as f:
				f.write(pretty_xml_str)
			print('The modified XMF files is written as:  ' + output_name)
	else:
		print(field_search + " element not found.")