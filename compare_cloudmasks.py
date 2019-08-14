#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''

Program reads cloudmasks either from different sensors or for different algorithms and compares them with each other

'''


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
import calendar
from sklearn.metrics import confusion_matrix

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
	ds = None
	return data, gt, srs_wkt

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
	
	
def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    day = 0
    
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
	print "jd ", jd
	date = "%s%02d%02d" % (y,month,jd)
	print "date julian", date
	
	return date
    
    
def getSensor(_filen):
	
	sensorList = ["MODIS","VIIRS","AATSR","SLSTR"]
	
	for _sens in sensorList:
		if _sens in _filen:
			_fsensor = _sens
	
	if "MOD35" in _filen:
		_fsensor = "MODIS"
	elif "MODIS" in _filen:
		_fsensor = "MODIS"
	elif "SLSTR" in _filen:
		_fsensor = "SLSTR"
	elif "AVHRR" in _filen:
		_fsensor = "AVHRR"
	else:
		print "Unknown sensor"
		sys.exit(1)
		
	regionList = ["Asia", "Europe", "NorthAmerica"]
	
	for _region in regionList:
		if _region in _filen:
			_fregion = _region
			
	if ("SCDA2.0" in _filen) or ("SCDA2_0" in _filen) or ("SCDA2.tif" in _filen) or ("SCDA2.mdi" in _filen):
		_alg = "SCDA2"
	elif ("SCDA2_1" in _filen) or ("SCDA2.1" in _filen):
		_alg = "SCDA2_1"
	elif "Cloud_CCI" in _filen:
		_alg = "Cloud_CCI"
	elif "METOP" in _filen:
		_alg = "metop"
	elif "NOAA-18" in _filen:
		_alg = "noaa18"
	elif "NOAA-19" in _filen:
		_alg = "noaa19"
	elif "NOAA-15" in _filen:
		_alg = "noaa15"
	elif ("NOAA-16" in _filen) or ('noaa16' in _filen):
		_alg = "noaa16"
	elif ("NOAA-17" in _filen) or ('noaa17' in _filen):
		_alg = "noaa17"
	elif "Cloud_CCI" in _filen:
		_alg = "Cloud_CCI"
	else :
		_alg = "product"
			
	return _fsensor, _fregion, _alg


def importBandForSensor(_sensor,_sfile,_alg):
	if _sensor == "MODIS" and _alg == "Cloud_CCI":
		_cloud = "%s" % _sfile
		#_time = _sfile.split(".")[2]
	elif _sensor == "MODIS" and _alg != "Cloud_CCI":
		_cloud = "%s,1" % _sfile
		#_time = _sfile.split(".")[2]
		
	elif _sensor == "MODIS" and _alg == "Cloud_CCI":
		_cloud = "%s" % _sfile
		#_time = _sfile.split(".")[2]

	elif _sensor == "VIIRS":
		_cloud = "%s,1" % _sfile
		
	elif _sensor == "SLSTR" and _alg != "Cloud_CCI":
		#_cloud = "%s,1" % _sfile 
		_cloud = "%s" % _sfile 
	elif _sensor == "SLSTR" and _alg == "Cloud_CCI":
		_cloud = "%s" % _sfile 
		
	elif _sensor == "AVHRR":
		_cloud = "%s" % _sfile
	
	elif _sensor == "AATSR":
		_cloud = "%s" % _sfile
		
	else:
		print "No such sensor"
		sys.exit(1)
		
	filelist = {"2003/Asia" : "20030301", "2003/Europe" : "20030303","2003/NorthAmerica" : "20030306", "2017/Asia" : "20170313", "2017/Europe" : "20170312", "2017/NorthAmerica" : "20170228"}
	for key in filelist:
		if key in _sfile:
			_date = filelist[key]
	
	return _cloud, _date

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":

	#  populating command line parser with options
	parser = optparse.OptionParser()
	
	parser.add_option("--f1",
				action="store", type="string", dest="file1", default = None,
				help="Specify the absolute path of the first file")
	
	parser.add_option("--f2",
				action="store", type="string", dest="file2", default = None,
				help="Specify the absolute path of the second file")
	
	parser.add_option("--extent",
			action="store", type="string", dest="extent", default = None,
			help="Specify the extent (xmin,ymin,xmax,ymax")
	
	parser.add_option("-o",
			action="store", type="string", dest="out_dir", default = None,
			help="Specify output directory")
	
	#Abfragen der command line options
	(options, args) = parser.parse_args()
	
	f1 = options.file1
	f2 = options.file2
	outdir = options.out_dir
	ext_import = options.extent
	ext = ext_import.replace(","," ")
	print ext
	
	sensor1,region1,alg1 = getSensor(f1)
	sensor2,region2,alg2 = getSensor(f2)

	if os.path.exists(outdir) != True:
		os.mkdir(outdir)
	os.chdir(outdir)
	
	cloud1,date1 = importBandForSensor(sensor1,f1,alg1)
	print "sensor 1 = ", sensor1, "algorithm = ", alg1
	cloud2,date2 = importBandForSensor(sensor2,f2,alg2)
	print "sensor 2 = ", sensor2, "algorithm = ", alg2
	
	out_f1 = "cloud_f1_%s_%s_%s.tif" % (sensor1,date1,region1)
	out_f2 = "cloud_f2_%s_%s_%s.tif" % (sensor2,date2,region2)
	
	if os.path.exists(out_f1) != True:
		cmd = "gdalwarp -tr 0.01 0.01 -te %s %s %s" % (ext, cloud1,out_f1)
		cmdexec(cmd)
	if os.path.exists(out_f2) != True:
		cmd = "gdalwarp -tr 0.01 0.01 -te %s %s %s" % (ext, cloud2,out_f2)
		cmdexec(cmd)
	
	#generate difference map for cloudmask
	if region1 == region2 and date1 == date2:
		out_diff = "cloud_diffmap_%s_%s_%s_%s_%s.tif" % (sensor1,sensor2,date1,date2,region1)
	else:
		out_diff = "cloud_diffmap_%s_%s_%s_%s_%s_%s.tif" % (sensor1,sensor2,date1,date2,region1,region2)
	out_f_product1 = "cloudmask_%s_%s_%s_product.tif" % (sensor1, region1,alg1)
	out_f_product2 = "cloudmask_%s_%s_%s_product.tif" % (sensor2, region2,alg2)
	out_f_cci1 = "cloudmask_f1_%s_%s_%s.tif" % (sensor1, region1,alg1)
	out_f_cci2= "cloudmask_f2_%s_%s_%s.tif" % (sensor2, region2,alg2)
	
	print out_f_cci1, out_f_cci2
	f1_read,gt3,srs_wkt3 =read_geo(out_f1,1)
	f2_read,gt4,srs_wkt4 = read_geo(out_f2,1)
	

	if sensor1 == "SLSTR" and alg1 == "product" and sensor2 == "MODIS" and alg2 == "product":
		f1_new = copy.deepcopy(f1_read)
		
		f1_new[numpy.where(f1_read < 3000)]=30
		f1_new[numpy.where(f1_read > 3000 )]=255
		f1_new[numpy.where(f1_read == 65535)] = 255
		f1_new[numpy.where(f1_read == 0)] = 255
		
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 100)] = 30
		f2_new[numpy.where(f2_read != 100)] = 255
		layer = f2_new - f1_read
	
		layer = f2_new - f1_new
		
	elif sensor1 == "SLSTR" and alg1 == "Cloud_CCI":
		print "SLSTR and CLoud_CCI"
		f1_new = copy.deepcopy(f1_read)
		f1_new[numpy.where(f1_read == 1)] = 30
		f1_new[numpy.where(f1_read != 1)] = 255
		layer = f2_read - f1_new
		
	elif sensor2 == "SLSTR" and alg2 == "Cloud_CCI":
		print "SLSTR and CLoud_CCI"
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 1)] = 30
		f2_new[numpy.where(f2_read != 1)] = 255
		layer = f2_new - f1_read
		
	elif sensor2 == "SLSTR" and alg2 == "product":
		f2_new = copy.deepcopy(f2_read)
		
		f2_new[numpy.where(f2_read < 3000)]=30
		f2_new[numpy.where(f2_read > 3000 )]=255
		f2_new[numpy.where(f2_read == 65535)] = 255
		f2_new[numpy.where(f2_read == 0)] = 255
	
		layer = f2_new - f1_read
		
	elif sensor2 == "MODIS" and alg2 == "Cloud_CCI":
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 1)] = 30
		f2_new[numpy.where(f2_read != 1)] = 255
		layer = f2_new - f1_read
	
		layer = f2_new - f1_read
		
	elif sensor2 == "MODIS" and alg2 == "product":
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 100)] = 30
		f2_new[numpy.where(f2_read != 100)] = 255
		layer = f2_new - f1_read
		
	elif sensor1 == "AVHRR" and alg1 == "Cloud_CCI" and sensor2 == "AVHRR" and alg2 =="noaa16":
        
 		f1_new = copy.deepcopy(f1_read)
		f1_new[numpy.where(f1_read == 1)] = 30
		f1_new[numpy.where(f1_read != 1)] = 255
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 1)] = 30
		f2_new[numpy.where(f2_read != 1)] = 255
		
		layer = f2_new - f1_new
		
	elif sensor1 == "AVHRR" and alg1 == "Cloud_CCI" and sensor2 == "AVHRR" and alg2 =="noaa17":
        
 		f1_new = copy.deepcopy(f1_read)
		f1_new[numpy.where(f1_read == 1)] = 30
		f1_new[numpy.where(f1_read != 1)] = 255
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 1)] = 30
		f2_new[numpy.where(f2_read != 1)] = 255
		
		layer = f2_new - f1_new

	elif (sensor2 == "AVHRR" and alg2 == "metop") or (sensor2 == "AVHRR" and alg2 == "noaa18") or (sensor2 == "AVHRR" and alg2 == "noaa19"):
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where((f2_read == 2) ^ (f2_read == 3))] = 30
		f2_new[numpy.where((f2_read == 0) ^ (f2_read == 1))] = 255
		layer = f2_new - f1_read
		
	elif ((sensor2 == "AVHRR" and alg2 == "noaa15") or (sensor2 == "AVHRR" and alg2 == "noaa16") or (sensor2 == "AVHRR" and alg2 == "noaa17")) and ((sensor1 == "AVHRR" and alg1 == "noaa15") or (sensor1 == "AVHRR" and alg1 == "noaa16") or (sensor1 == "AVHRR" and alg1 == "noaa17")):
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where((f2_read == 2) ^ (f2_read == 3))] = 30
		f2_new[numpy.where((f2_read == 0) ^ (f2_read == 1))] = 255
		
		f1_new = copy.deepcopy(f1_read)
		f1_new[numpy.where((f1_read == 2) ^ (f1_read == 3))] = 30
		f1_new[numpy.where((f1_read == 0) ^ (f1_read == 1))] = 255
		layer = f2_new - f1_new
	elif sensor1 == "AVHRR" and alg1 == "Cloud_CCI" and sensor2 == "AVHRR" and alg2 =="Cloud_CCI":
        
 		f1_new = copy.deepcopy(f1_read)
		f1_new[numpy.where(f1_read == 1)] = 30
		f1_new[numpy.where(f1_read != 1)] = 255
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 1)] = 30
		f2_new[numpy.where(f2_read != 1)] = 255
		layer = f2_new - f1_new		
	elif sensor1 == "AVHRR" and alg1 == "Cloud_CCI":
		f2_new = copy.deepcopy(f2_read)
		f2_new[numpy.where(f2_read == 1)] = 30
		f2_new[numpy.where(f2_read != 1)] = 255
		layer = f2_new - f1_read

		
	else:
		layer = f2_read - f1_read
		
	layer[numpy.where(layer == 0)]= 0 # the same values
	layer[numpy.where(layer != 0)] = 1
	
	
	if os.path.exists(out_diff) != True:
		make_geo(layer,out_diff,gt3,srs_wkt3,1)
	
	if alg1 == "product":
		make_geo(f1_new,out_f_product1,gt4,srs_wkt4,1)
	if alg2 == "product":
		make_geo(f2_new,out_f_product2,gt4,srs_wkt4,1)
	if alg1 == "Cloud_CCI":
		make_geo(f1_new,out_f_cci1,gt4,srs_wkt4,1)
	if alg2 == "Cloud_CCI":
		make_geo(f2_new,out_f_cci2,gt4,srs_wkt4,1)
	if alg2 == "metop" or alg2 == "noaa18" or alg2 == "noaa19":
		make_geo(f2_new,out_f_cci2,gt4,srs_wkt4,1)
	if alg2 == "noaa15" or alg2 == "noaa16" or alg2 == "noaa17":
		make_geo(f2_new,out_f_cci2,gt4,srs_wkt4,1)
	if alg1 == "noaa15" or alg1 == "noaa16" or alg1 == "noaa17":
		make_geo(f1_new,out_f_cci1,gt4,srs_wkt4,1)
		
	#confusion matrix
	if os.path.exists(out_f_product1):
		f1_stat = numpy.array(f1_new).flatten()
	elif os.path.exists(out_f_cci1): 
		f1_stat = numpy.array(f1_new).flatten()
	else:
		f1_stat = numpy.array(f1_read).flatten()
	
	if os.path.exists(out_f_product2):
		f2_stat = numpy.array(f2_new).flatten()
	elif os.path.exists(out_f_cci2): 
		f2_stat = numpy.array(f2_new).flatten()
	else:
		f2_stat = numpy.array(f2_read).flatten()
		
	x = confusion_matrix(f1_stat,f2_stat)
	print x
	
	true_p = 100*float(x[0][0])/float((x[0][0]+x[0][1]))
	false_p = 100*float(x[0][1])/float((x[0][0]+x[0][1]))
	false_n = 100*float(x[1][0])/float((x[1][0]+x[1][1]))
	true_n = 100*float(x[1][1])/float((x[1][0]+x[1][1]))
	print "%.2f" % true_p, "%.2f" % false_p, "%.2f" % false_n, "%.2f" % true_n
	
	ov_acc = (true_p + true_n)/(true_p + false_p + false_n + true_n)
	
	logfile = "confusionmatrix_%s_%s_%s_%s.csv" % (sensor1,alg1,sensor2,alg2)
	if os.path.exists(logfile):
		cmdexec("rm %s" % logfile)
	
	f=open(logfile,"w")
	f.write("%.2f %.2f" % (true_p,false_p))
	f.write(",")
	f.write("%.2f %.2f" % (false_n,true_n))
	f.write("\n")
	f.write("%.2f" % ov_acc)
	f.close()
