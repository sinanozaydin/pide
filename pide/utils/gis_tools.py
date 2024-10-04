import math
from pyproj import Proj, transform
import numpy as np

def get_utm_zone_number(longitude):

	"""
	A function to get utm zone number
	Input:
	float: longitude
	
	Output:
	UTM Zone Number
	"""

	zone_number = math.ceil((longitude + 180.0) / 6.0)
	
	return zone_number

def lat_lon_to_utm(latitude, longitude, zone_number = None):
	"""
	A function to convert latitude longitude coordinates in WGS84 projection to UTM coordinates.
	
	Input:
	float or array: latitude
	float or array: longitude
	
	Output:
	UTM Coordinates
	"""
	wgs84 = Proj(proj='latlong', datum='WGS84')
	
	if zone_number == None:
		
		if isinstance(longitude, (float,int,np.float64,np.int64)):
			zone_number = get_utm_zone_number(longitude=longitude)
		else:
			raise ValueError('You have to enter a valid zone number! You can get this by using get_utm_zone_number method')
	
	utm = Proj(proj='utm', zone=zone_number, datum='WGS84')

	utm_x, utm_y = transform(wgs84, utm, longitude, latitude)
	
	return utm_x, utm_y
	
def utm_to_lat_lon(utm_x, utm_y, zone_number):
	"""
	A function to convert UTM coordinates to latitude longitudes in WGS84 projectin.
	"""
	utm = Proj(proj='utm', zone=zone_number, datum='WGS84', ellps='WGS84')

	longitude, latitude = transform(utm, Proj(proj='latlong', datum='WGS84'), utm_x, utm_y)
	
	return latitude, longitude