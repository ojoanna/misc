#!/usr/bin/env python

# -*- coding: utf-8 -*-

import time, os, optparse, sys
import os
import numpy, pylab
import osgeo.gdal as gdal
import osgeo.gdalconst as gdalconst
import osr
import pdb
import copy
import glob
import gdalnumeric
import datetime
import ogr

def cmdexec(cmd):
	print cmd
	try:
		#retcode = call(cmd, shell=True)
		retcode = os.system(cmd)
		#retcode = 0
		if retcode != 0:
			print >>sys.stderr, "Child was terminated by signal", -retcode
			#sys.exit(1)
		#else:
			#print >>sys.stderr, "Child returned", retcode
	except OSError, e:
		print >>sys.stderr, "Execution failed:\n<%s>", e, cmd
		sys.exit(1)


def read_geo(name, layer):
	ds=None
	ds = gdal.Open(name, gdal.GA_ReadOnly)
	data = ds.GetRasterBand(layer).ReadAsArray()
	gt = ds.GetGeoTransform()
	srs_wkt = ds.GetProjection()
	pixel_size=ds.GetGeoTransform()[1]
	ds = None
	return data, gt, srs_wkt, pixel_size

def get_projection(name):
	ds=None
	ds=gdal.Open(name, gdal.GA_ReadOnly)
	prj=ds.GetProjection()
	srs=osr.SpatialReference(wkt=prj)
	utm=srs.GetUTMZone()
	
	return utm

def make_geo(data, output_name, gt, srs_wkt, layer):
	#make Geotiff of difference and quiklook png
	#rm old file_
	try:
		os.system('rm ' + output_name)
	except:
		'No such file or directory'
		pass

	ds_save = None
	Ny, Nx = data.shape
	driver = gdal.GetDriverByName('GTiff')
	ds_save = driver.Create(output_name, Nx, Ny, int(layer),  gdal.GDT_Float32)#write range to tif
	ds_save.SetGeoTransform(gt)
	ds_save.SetProjection(srs_wkt)
	ds_save.GetRasterBand(1).WriteArray(data)
	ds_save.GetRasterBand(1).SetNoDataValue(255)
	ds_save = None
	
def get_coordinates(_shp):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(_shp, 0)
	layer = dataSource.GetLayer()
	
	for feature in layer:
		_geom = feature.GetGeometryRef()
		
	_envelope = _geom.GetEnvelope()
	_xmin=_envelope[0]
	_xmax=_envelope[1]
	_ymin=_envelope[2]
	_ymax=_envelope[3]
	
	print _xmin,_xmax,_ymin,_ymax
	_coordstr = "%s %s %s %s" % (str(_xmin)[:4],str(_ymin)[:4],str(_xmax)[:4],str(_ymax)[:4])
	
	return _coordstr
	
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":

	#  populating command line parser with options
	parser = optparse.OptionParser()
	
	parser.add_option("--date",
				   action="store", type="string", dest="date", default = None,
				   help="Specify the date, format YYYYMMDD")
		
	parser.add_option("--SAR",
				   action="store", type="string", dest="sar_im_dir", default = None,
				   help="Specify wet snowmap directory created with SAR data.")
	
	parser.add_option("--OPT",
				   action="store", type="string", dest="opt_im_dir", default = None,
				   help="Specify dry&wet snowmap directory created with optical data")
	parser.add_option("--output_dir",
				   action="store", type="string", dest="out_dir", default = None,
				   help="Specify output directory.")
	parser.add_option("--track",
				   action="store", type="string", dest="track", default = None,
				   help="Specify track number.")
	#parser.add_option("--utm",
				   #action="store", type="string", dest="utm", default = None,
				   #help="Specify utm.")
	

	#Abfragen der command line options
	(options, args) = parser.parse_args()

	options, args = parser.parse_args(sys.argv[1:])
	date = options.date
	sarImageDir = options.sar_im_dir
	optImageDir = options.opt_im_dir
	outDir = options.out_dir
	track = options.track
	#coordstr = options.cord
	#utm_n=options.utm
	
	kmlDir = "/mnt/h5/projects/S14SCI_snow/NtoS_track168/kml"
	shapefile = "%s/AOI_Europe_track%s.shp" % (kmlDir,track)
	coord = get_coordinates(shapefile)

	
	
	#set parameters 
	
	year = date[:4]
	monthday = date[4:]
	dateInt=[int(date[:4]), int(date[4:6]),int(date[6:])]
	yddate = "%0*d" % (3,(datetime.datetime(dateInt[0],dateInt[1],dateInt[2]).timetuple().tm_yday))
		
	#sarImageRaw = "%s/%s/m16x4/old/sm_avq_%s.img" % (sarImageDir,date,date)
	#if os.path.exists(sarImageRaw) != True:
	sarImageRaw = "%s/%s/m16x4/sm_avq_%s.img" % (sarImageDir,date,date)
	print sarImageRaw
	optdatedir = "%s/%s/%s/LatLon_FSCA*.tif" % (optImageDir, year, monthday)  #for Johanna's assimilated product
	optImageRaw = glob.glob(optdatedir)[0]
	
	output_geo="wetSnowMap100m_%s.tif" % date
	#output_geo2="FSC_wetDrySnowMap_OPT_%s.tif" % date

	#coordlist=coordstr.split(",")
	#coord=(" ").join(coordlist)
	#coord="-11.0 35.0 50.0 72.0"
	
	if os.path.exists(sarImageRaw) == True:

		#optdatedir = "%s/%s/%s/FSC_0.005deg_*.tif" % (optImageDir, year, yddate)  #for cryoland in SNOW_PROCESSING directory products

		try:
			os.mkdir(outDir)
		except:
			pass
		
		outDirDate = "%s/%s/" % (outDir, date)
		if os.path.exists(outDirDate) != True:
			cmdexec("mkdir %s" % outDirDate)
		
		os.chdir(outDirDate)
		utm_n=get_projection(sarImageRaw)
		
		optImage2 = os.path.basename(optImageRaw)[:-4]+"_cut2.tif"
		optImage = os.path.basename(optImageRaw)[:-4]+"_cut.tif"
		if os.path.isfile(optImage2) == True:
			cmdexec("rm %s" % optImage2)
		cmd = "gdalwarp -te %s %s %s" % (coord,optImageRaw, optImage2)
		cmdexec(cmd)		
		
		if os.path.isfile(optImage) == True:
			cmdexec("rm %s" % optImage)
		cmd = "gdalwarp -tr 0.001 0.001 %s %s" % (optImage2,optImage)
		cmdexec(cmd)
		
		sarImage = os.path.basename(sarImageRaw)[:-4]+"_4326.tif"
		sarImage2 = os.path.basename(sarImageRaw)[:-4]+"_4326_2.tif"

		if os.path.isfile(sarImage2) == True:
			cmdexec("rm %s" % sarImage2)
		cmd = "gdalwarp -s_srs EPSG:326%s -t_srs EPSG:4326 -tr 0.001 0.001 -te %s %s %s" % (utm_n, coord, sarImageRaw, sarImage2)
		cmdexec(cmd)
		
		if os.path.isfile(sarImage) == True:
			cmdexec("rm %s" % sarImage)
		cmd = "gdalwarp -tr 0.001 0.001 -te %s %s %s" % (coord, sarImage2, sarImage)
		cmdexec(cmd)
		
		
		optImageData, gt1, proj1,pix1 = read_geo(optImage,1)

		#remove wet snow where 
		sarImageData,gt2, proj2,pix2 = read_geo(sarImage,1)
		
		layer = copy.deepcopy(sarImageData)
		#layer[numpy.where(sarImageData==216)]=sarImageFSC[numpy.where(sarImageData==216)]
		#layer[numpy.where(sarImageData==209)]=100
		#layer[numpy.where((optImageData<170) & (optImageData>=100) & (sarImageData==209))]=100
		layer[numpy.where((optImageData<170) & (optImageData>=100) & (sarImageData==216))]=209

		
		optImageData[numpy.where(layer==255)]=255
		optImageData[numpy.where(optImageData==30)]=255
		
		make_geo(layer,output_geo,gt1,proj1,1)
		#make_geo(optImageData,output_geo2,gt1,proj1,1)
		
		#remove temporary files
		cmd = "rm %s %s %s %s" % (sarImage, optImage, sarImage2, optImage2)
		cmdexec(cmd)
		
		print "GeoTiff created!"
	else:
		print "Wet snow map is not available for this date!"
		pass
	
	###count wetsnow pixels
	#raster_file = gdalnumeric.LoadFile(output_geo)
	#pixel_count = (raster_file == 216).sum()
	#print "Wet snow pixels = ", pixel_count
	
