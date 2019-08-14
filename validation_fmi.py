#!/usr/bin/env python

'''
Program validates wet snow product created by our company comparing to the files that we got from Finland.
'''

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

def make_geo(data, output_name, gt, srs_wkt, layer):
	#make Geotiff of difference and quiklook png
	#rm old file_
	try:
		os.system('rm ' + output_name)
	except:
		'No such file or directory'
		pass

	ds_save = None
	#pdb.set_trace()
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
	_coordstr = "%s %s %s %s" % (_xmin,_ymin,_xmax,_ymax)
	
	return _coordstr



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":

	#  populating command line parser with options
	parser = optparse.OptionParser()
	
	parser.add_option("--date",
				   action="store", type="string", dest="date", default = None,
				   help="Specify the date, format YYYYMMDD")
		
	parser.add_option("--enveo",
				   action="store", type="string", dest="enveo_dir", default = None,
				   help="Specify wet snowmap directory created with SAR data.")
	
	parser.add_option("--fmi",
				   action="store", type="string", dest="fmi_dir", default = None,
				   help="Specify dry&wet snowmap directory created with optical data")
	
	parser.add_option("--track",
				   action="store", type="string", dest="track", default = None,
				   help="Specify the track, format XXX")
				   
				   
				   
	#Abfragen der command line options
	(options, args) = parser.parse_args()

	options, args = parser.parse_args(sys.argv[1:])
	date = options.date
	enveoImageDir = options.enveo_dir
	fmiImageDir = options.fmi_dir
	track=options.track
	
	if os.path.exists(enveoImageDir)!=True:
		print "wrong directory!"
		sys.exit(1)
		
	if os.path.exists(fmiImageDir)!=True:
		print "wrong directory!"
		sys.exit(1)
		
	root_dir = "/mnt/h1/projects/S14SCI_snow/data_fmi/processing"
	trackdir="%s/track%s/%s" % (enveoImageDir, track, date)
	ws_file_raw = "%s/FSC_wetSnowMap_%s.tif" % (trackdir,date)
	ws_file_temp = "FSC_wetSnowMap_%s_res_track%s.tif" % (date,track)
	fmi_file_raw = "%s/%s_day_HSAF_SS.tif" % (fmiImageDir, date)
	fmi_file_temp = "%s_day_HSAF_SS_res_track%s.tif" % (date,track)
	
	os.chdir(root_dir)
		
	if os.path.exists(trackdir) != True:
		print "Wet snow data is not available for this date!"
		sys.exit(1)
	
	shapefile = "/mnt/h5/projects/S14SCI_snow/NtoS_track168/kml/AOI_Europe_track%s.shp" % track
	coord = get_coordinates(shapefile)
	
	if os.path.exists(fmi_file_temp) != True:
		cmd = "gdalwarp -tr 0.25 0.25 -te %s %s %s" % (coord, fmi_file_raw, fmi_file_temp)
		cmdexec(cmd)
	
	if os.path.exists(ws_file_temp) != True:
		cmd ="gdalwarp -tr 0.25 0.25 -te %s %s %s" % (coord, ws_file_raw, ws_file_temp)
		cmdexec(cmd)
	
	ws_data,gt1,proj1,pixel = read_geo(ws_file_temp,1)
	fmi_data,gt2,proj2,pixel2 = read_geo(fmi_file_temp,1)
	
	ws_data[numpy.where((ws_data>100) & (ws_data<201))]=200
	ws_data[numpy.where((ws_data==200) & (fmi_data==42))]=255
	fmi_data[numpy.where(ws_data==255)]=255
	layer = copy.deepcopy(ws_data)
	layer[numpy.where((fmi_data==25) & (ws_data==200))]=99
	layer[numpy.where((fmi_data!=25) & (ws_data!=200))]=255
	
	output_fmi_data = "fmi_data_200_track%s_%s.tif" % (track,date)
	make_geo(fmi_data,output_fmi_data,gt1,proj1,1)
	
	output_ws_data = "ws_data_200_track%s_%s.tif" % (track,date)
	make_geo(ws_data,output_ws_data,gt1,proj1,1)
	
	output_geo ="layer_track%s_%s.tif" % (track,date)
	make_geo(layer,output_geo,gt1,proj1,1)
		  
	##count wetsnow pixels
	raster_file = gdalnumeric.LoadFile(output_geo)
	raster_file_enveo = gdalnumeric.LoadFile(output_ws_data)
	raster_file_fmi = gdalnumeric.LoadFile(output_fmi_data)
	pixel_count = (raster_file == 99).sum()
	pixel_count_enveo = (raster_file_enveo == 200).sum()
	pixel_count_fmi = (raster_file_fmi == 25).sum()
	print "Wet snow pixels common= ", pixel_count
	print "Wet snow pixels enveo= ", pixel_count_enveo
	print "Wet snow pixels fmi= ", pixel_count_fmi
	
	#write the results to the logfile
	logfile="summary.csv"
	str_data=[track,date,str(pixel_count),str(pixel_count_enveo),str(pixel_count_fmi)]
	str_csv=','.join(str_data)
	if os.path.exists(logfile)!=True:
		fl=open(logfile,'w')
		fl.write("track,date,common,enveo,fmi")
		fl.write("\n")
		fl.write(str_csv)
		fl.close()
	else:
		fl=open(logfile,'a')
		fl.write("\n")
		fl.write(str_csv)
		fl.close()
	cmdexec("rm %s %s" % (ws_file_temp,fmi_file_temp))
	
	
