import os
import argparse
import pdb
from osgeo import ogr

'''

Read the shapefile's coordinates and convert them into the way that can be used on https://scihub.copernicus.eu/ or https://search.asf.alaska.edu/

'''

def read_geometry(_shapefile):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(_shapefile,0)
	layer = dataSource.GetLayer()
	for feature in layer:
		geom = feature.GetGeometryRef()
		
	_polygon = str(geom)
	
	_coordinates = _polygon.strip("POLYGON ((").strip("))")

	
	return _coordinates


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=("Read the shapefile's coordinates and convert them into the way that can be used on https://scihub.copernicus.eu/ or https://search.asf.alaska.edu/"))
        
	parser.add_argument('-f', '--shapefile', metavar='shapefile', action='store', type=str, required=False, dest='shapefile',
                                        help='specify shapefile of the track')
	parser.add_argument('-t', '--track', metavar='track', action='store', type=str, required=False, dest='track',
                                        help='specify number of the track')
	
	shp_dir = "/mnt/h5/projects/S14SCI_snow/NtoS_track168/kml/"
	args = parser.parse_args()
	
	if args.shapefile:
		shapefile = args.shapefile
	else:
		track = args.track
		shapefile = "%s/AOI_Europe_track%s.shp" % (shp_dir,track)
		print shapefile
	

	coords = read_geometry(shapefile)
	newcords2 = coords.replace(" 0","")
	print newcords2
	
	newll=coords.replace(' ',',')
	newlll=newll.strip("\n")
	#newll=line.strip("\n")
	file_list = newlll.split(",")
	new_list=[]
	#pdb.set_trace()
	for nnumber in file_list:
		if nnumber != "0":
			try:
				x = round(float(nnumber),2)
				new_list.append(x)
			except:
				continue
		
	print new_list[:-2]


