#!/usr/bin/env python


'''

The program reads the necessary bands depending on the sensor and creates NDSI and NDVI difference map as well as figures of their dependency.


'''
# -*- coding: utf-8 -*-

import time, os, optparse, sys
import numpy
import osgeo.gdal as gdal
import osgeo.gdalconst as gdalconst
import osr
import pdb
import copy
import glob
import gdalnumeric
import datetime
import calendar
import matplotlib.pyplot as plt
from scipy import stats

import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms

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
			
	regionList = ["Asia", "Europe", "NorthAmerica"]
	for _region in regionList:
		if _region in _filen:
			_fregion = _region
			
	return _fsensor, _fregion


def importBandForSensor(_sensor,_sfile):
	
	if _sensor=="MODIS":
		_r550 = "%s,_ev_500_refsb_b4" % _sfile
		_r1600 = "%s,_ev_500_refsb_b6" % _sfile
		_r650 = "%s,_ev_250_refsb_b1" % _sfile
		_r850 = "%s,_ev_250_refsb_b2" % _sfile
		#_yyyy = _sfile.split(".")[1][1:5]
		#_ddd = _sfile.split(".")[1][5:8]
		#pdb.set_trace()
		#if _ddd[0] == "0":
			#_ddd = _ddd[1:]
		#_date = JulianDate_to_MMDDYYY(int(_yyyy),int(_ddd))
		#print "date MODIS = ", _date
		
	elif _sensor == "VIIRS":
		_r550 = "%s,_m04" % _sfile
		_r650 = "%s,_m05" % _sfile
		_r850 = "%s,_m07" % _sfile
		_r1600 = "%s,_m10" % _sfile
		
		#_date = _sfile.split("_")[1][1:9]
		
	elif _sensor == "SLSTR":
		_r550 = "%s,_toa_reflectance_s1_an" % _sfile 
		_r650 = "%s,_toa_reflectance_s2_an" % _sfile 
		_r850 = "%s,_toa_reflectance_s3_an" % _sfile
		_r1600 = "%s,_toa_reflectance_s5_an" % _sfile
		#_date = _sfile.split("____")[1][:8]
		
		
	elif _sensor == "AATSR":
		_r550 = "%s" % _sfile
		_r650 = "%s" % _sfile
		_r850 = "%s" % _sfile
		_r1600 = "%s" % _sfile
		#_date = _sfile.split("/")[-1].split("_")[5][6:14]
		#print "aatsr date ", _date
	else:
		print "No such sensor"
		sys.exit(1)
		
	filelist = {"2003/Asia" : "20030301", "2003/Europe" : "20030303","2003/NorthAmerica" : "20030306", "2017/Asia" : "20170313", "2017/Europe" : "20170312", "2017/NorthAmerica" : "20170228"}
	for key in filelist:
		if key in _sfile:
			_date = filelist[key]
			
	return _r550, _r650, _r850, _r1600, _date

def calculateStatistics(_f1,_f2,_statfile,_type):
	_f1_a = numpy.array(_f1).flatten()
	_f2_a = numpy.array(_f2).flatten()
	
	_cor = numpy.corrcoef(_f1_a,_f2_a)[0][1]
	print "correlation coefficient = ",_cor
	#print f2_a
	
	_bias = (_f2_a-_f1_a).sum()/float(_f1_a.size)
	print "bias = ", _bias
	
	_RMSD = numpy.sqrt( ((_f2_a-_f1_a)*(_f2_a-_f1_a)).sum() / float(_f1_a.size) )
	print "RMSD = ", _RMSD
	
	#_a=numpy.array(_f1_a,_f2_a)
	#_stdv = numpy.std(_f1_a,_f2_a)
	#print "standar deviation = ",_stdv
	
	f=open(_statfile,"a")
	f.write(_type)
	f.write("\n")
	f.write("correlation coefficient = %s" % _cor)
	f.write("\n")
	f.write("bias = %s" % _bias)
	f.write("\n")
	f.write("RMSD = %s" % _RMSD)
	f.write("\n")
	#f.write("standar deviation = " % _stdv)
	#f.write("\n")
	f.write("\n")
	f.close()
	
	
	return _cor,_bias,_RMSD
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
	
	parser.add_option("--id",
			action="store", type="string", dest="idname", default = None,
			help="Specify special id for the file")
	
	
	#Abfragen der command line options
	(options, args) = parser.parse_args()
	
	f1 = options.file1
	f2 = options.file2
	outdir = options.out_dir
	ext_import = options.extent
	ext = ext_import.replace(","," ")
	idname = options.idname
	print "AOI = ", ext
	
	
	llc_x = "E%s" % ext_import.split(",")[0].replace(".","")
	if "-" in llc_x:
		llc_x = "W%s" % llc_x[2:]

	llc_y = "N%s" %  ext_import.split(",")[1].replace(".","")
	if "-" in llc_x:
		llc_y = "S%s" % llc_y[2:]

	urc_x = "E%s" % ext_import.split(",")[2].replace(".","")
	if "-" in urc_x:
		urc_x = "W%s" % urc_x[2:]

	urc_y = "N%s" %  ext_import.split(",")[3].replace(".","")
	if "-" in urc_y:
		urc_y = "S%s" % urc_y[2:]
	
	sensor1,region1 = getSensor(f1)
	sensor2,region2 = getSensor(f2)
	
	#import band for specific sensors
	f1R550mdi, f1R650mdi, f1R850mdi, f1R1600mdi,date1 = importBandForSensor(sensor1,f1)
	f2R550mdi, f2R650mdi, f2R850mdi, f2R1600mdi,date2 = importBandForSensor(sensor2,f2)	
	print "sensor 1 = ", sensor1
	print "sensor 2 = ", sensor2
	
	#create output directory
	out_region = "%s/%s/%s/%s_%s_%s_%s_%s" % (outdir,date1[:4],region1,idname,llc_y,urc_y,llc_x,urc_x)
	if os.path.exists(out_region) != True:
		os.makedirs(out_region)
		
	os.chdir(out_region)
	
	#extract the bands from the mdi files and cut them to the desired extent to reduce the time of reading the data later
	f1R550 = "%s_f1r550.tif" % os.path.basename(f1R550mdi).split(",")[0][:-4]
	f1R650 = "%s_f1r650.tif" % os.path.basename(f1R650mdi).split(",")[0][:-4]
	f1R850 = "%s_f1r850.tif" % os.path.basename(f1R850mdi).split(",")[0][:-4]
	f1R1600 = "%s_f1r1600.tif" % os.path.basename(f1R1600mdi).split(",")[0][:-4]
	f2R550 = "%s_f2r550.tif" % os.path.basename(f2R550mdi).split(",")[0][:-4]
	f2R650 = "%s_f2r650.tif" % os.path.basename(f2R650mdi).split(",")[0][:-4]
	f2R850 = "%s_f2r850.tif" % os.path.basename(f2R850mdi).split(",")[0][:-4]
	f2R1600 = "%s_f2r1600.tif" % os.path.basename(f2R550mdi).split(",")[0][:-4]
	
	if os.path.exists(f1R550) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f1R550mdi,f1R550)
		cmdexec(cmd)
	if os.path.exists(f1R650) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f1R650mdi,f1R650)
		cmdexec(cmd)
	if os.path.exists(f1R850) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f1R850mdi,f1R850)
		cmdexec(cmd)
	if os.path.exists(f1R1600) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f1R1600mdi,f1R1600)
		cmdexec(cmd)
	if os.path.exists(f2R550) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f2R550mdi,f2R550)
		cmdexec(cmd)
	if os.path.exists(f2R650) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f2R650mdi,f2R650)
		cmdexec(cmd)
	if os.path.exists(f2R850) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f2R850mdi,f2R850)
		cmdexec(cmd)
	if os.path.exists(f2R1600) != True:
		cmd = "gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -te %s %s %s" % (ext, f2R1600mdi,f2R1600)
		cmdexec(cmd)

	#read data
	if sensor1 == "AATSR":
		f1R550data, gt1, srs_wkt1= read_geo(f1R550, 7)
		f1R650data, gt1, srs_wkt1 = read_geo(f1R650, 6)
		f1R850data, gt1, srs_wkt1 = read_geo(f1R850, 5)
		f1R1600data, gt1, srs_wkt1= read_geo(f1R1600, 4)
		
		f2R550data, gt2, srs_wkt2= read_geo(f2R550, 1)
		f2R650data, gt2, srs_wkt2= read_geo(f2R650, 1)
		f2R850data, gt2, srs_wkt2= read_geo(f2R850, 1)
		f2R1600data, gt2, srs_wkt2= read_geo(f2R1600, 1)
	elif sensor2 == "AATSR":
		f1R550data, gt1, srs_wkt1= read_geo(f1R550, 1)
		f1R650data, gt1, srs_wkt1 = read_geo(f1R650, 1)
		f1R850data, gt1, srs_wkt1 = read_geo(f1R850, 1)
		f1R1600data, gt1, srs_wkt1= read_geo(f1R1600, 1)
		
		f2R550data, gt2, srs_wkt2= read_geo(f2R550, 7)
		f2R650data, gt2, srs_wkt2 = read_geo(f2R650, 6)
		f2R850data, gt2, srs_wkt2 = read_geo(f2R850, 5)
		f2R1600data, gt2, srs_wkt2= read_geo(f2R1600, 4)
	else:
		f1R550data, gt1, srs_wkt1= read_geo(f1R550, 1)
		f1R650data, gt1, srs_wkt1 = read_geo(f1R650, 1)
		f1R850data, gt1, srs_wkt1 = read_geo(f1R850, 1)
		f1R1600data, gt1, srs_wkt1= read_geo(f1R1600, 1)
		
		f2R550data, gt2, srs_wkt2= read_geo(f2R550, 1)
		f2R650data, gt2, srs_wkt2= read_geo(f2R650, 1)
		f2R850data, gt2, srs_wkt2= read_geo(f2R850, 1)
		f2R1600data, gt2, srs_wkt2= read_geo(f2R1600, 1)
	
	#the band S5 of SLSTR has to be multiplied by the factor 1.12
	if sensor1 == "SLSTR":
		f1R1600data=1.12*f1R1600data
	elif sensor2 == "SLSTR":
		f2R1600data=1.12*f2R1600data

	f1R550data[numpy.where(f1R550data>1000)]=numpy.nan
	f1R650data[numpy.where(f1R650data>1000)]=numpy.nan
	f1R850data[numpy.where(f1R850data>1000)]=numpy.nan
	f1R1600data[numpy.where(f1R1600data>1000)]=numpy.nan
	f2R550data[numpy.where(f2R550data>1000)]=numpy.nan
	f2R1600data[numpy.where(f2R1600data>1000)]=numpy.nan
	
	#calculate NDSI and NDVI
	ndsi_f1 = (f1R550data-f1R1600data)/(f1R550data+f1R1600data)
	ndsi_f2 = (f2R550data-f2R1600data)/(f2R550data+f2R1600data)
	
	ndvi_f1 = (f1R850data - f1R650data)/(f1R850data + f1R650data)
	ndvi_f2 = (f2R850data - f2R650data)/(f2R850data + f1R650data)
	

	out_f1 = "ndsi2_f1_%s_%s_%s_%s%s%s%s.tif" % (sensor1,date1,region1,llc_y,urc_y,llc_x,urc_x)
	out_f2 = "ndsi2_f2_%s_%s_%s_%s%s%s%s.tif" % (sensor2,date2,region2,llc_y,urc_y,llc_x,urc_x)
	
	out_f1_ndvi = "ndvi2_f1_%s_%s_%s_%s%s%s%s.tif" % (sensor1,date1,region1,llc_y,urc_y,llc_x,urc_x)
	out_f2_ndvi = "ndvi2_f2_%s_%s_%s_%s%s%s%s.tif" % (sensor2,date2,region2,llc_y,urc_y,llc_x,urc_x)
	#out_f1_temp = "%s_temp.tif" % out_f1[:-4]
	#out_f2_temp = "%s_temp.tif" % out_f2[:-4]
	
	if os.path.exists(out_f1) != True:
		make_geo(ndsi_f1,out_f1,gt1,srs_wkt1,1)
	if os.path.exists(out_f2) != True:
		make_geo(ndsi_f2,out_f2,gt2,srs_wkt2,1)
	if os.path.exists(out_f1_ndvi) != True:
		make_geo(ndvi_f1,out_f1_ndvi,gt1,srs_wkt1,1)
	if os.path.exists(out_f2_ndvi) != True:
		make_geo(ndvi_f2,out_f2_ndvi,gt2,srs_wkt2,1)
	
	#generate difference map for NDSI
	if region1 == region2 and date1 == date2:
		out_diff = "ndsi_diffmap2_%s_%s_%s_%s_%s%s%s%s.tif" % (sensor1,sensor2,date1,region1,llc_y,urc_y,llc_x,urc_x)
	else:
		out_diff = "ndsi_diffmap2_%s_%s_%s_%s_%s_%s_%s%s%s%s.tif" % (sensor1,sensor2,date1,date2,region1,region2,llc_y,urc_y,llc_x,urc_x)
	

	layer = ndsi_f1 - ndsi_f2
	
	if os.path.exists(out_diff) != True:
		make_geo(layer,out_diff,gt1,srs_wkt1,1)
	


	# STATISTICS 
	#the extent and the files used 
	stat_file = "statistics_%s_%s_%s_%s.txt" % (sensor1,sensor2,region1,date1)
	
	ffstat = open(stat_file,"w")
	ffstat.write("Calculations for the files:\n")
	ffstat.write("file1 : %s\n" % f1)
	ffstat.write("file2 : %s\n\n" % f2)
	ffstat.write("extent (xmin,ymin,xmax,ymax) = %s\n" % (ext_import))
	ffstat.close()
	
	#plots
	fig = plt.figure( figsize=(7,7))
	fig.subplots_adjust(hspace=0.5, wspace=0.3)
	#NDSI 
	cor_ndsi,bias_ndsi,RMSD_ndsi = calculateStatistics(ndsi_f1,ndsi_f2,stat_file,"NDSI")
	ax1 = fig.add_subplot(2,2,1)
	#ax1.scatter( ndsi_f1, ndsi_f2, c="red",s=10,marker=".")
	ax1.set_title("NDSI")
	ax1.set_ylabel(sensor2)
	ax1.set_xlabel(sensor1)
	ax1.set_xlim(-0.5,1)
	ax1.set_ylim(-0.5,1)
	ax1.annotate("r = %.3f" % cor_ndsi,xy=(0.25, 0.94),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
	#regression line
	x_nsdi = numpy.array(ndsi_f1).flatten()
	y_nsdi = numpy.array(ndsi_f2).flatten()
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_nsdi,y_nsdi)
	line_nsdi = slope*x_nsdi+intercept
	ax1.plot(x_nsdi,y_nsdi,'ro', x_nsdi, line_nsdi,"b-",linewidth=0.5,markersize=1)

	
	#NDVI 
	cor_ndvi,bias_ndvi,RMSD_ndvi = calculateStatistics(ndvi_f1,ndvi_f2,stat_file,"NDVI")
	ax2 = fig.add_subplot(2,2,2)
	#ax2.scatter( ndvi_f1, ndvi_f2, c="green",s=10,marker=".")
	ax2.set_title("NDVI")
	ax2.set_xlim(-0.5,0.8)
	ax2.set_ylim(-0.5,0.8)
	ax2.set_xlabel(sensor1)
	ax2.annotate("r = %.3f" % cor_ndvi,xy=(0.25, 0.94),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
	
	#regression line
	x_ndvi = numpy.array(ndvi_f1).flatten()
	y_ndvi = numpy.array(ndvi_f2).flatten()
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_ndvi,y_ndvi)
	line_ndvi = slope*x_ndvi+intercept
	ax2.plot(x_ndvi,y_ndvi,'go', x_ndvi, line_ndvi,"b-",linewidth=0.5,markersize=1)


	#R550
	cor_r550, bias_r550,RMSD_r550 = calculateStatistics(f1R550data,f2R550data,stat_file,"R550")
	ax3 = fig.add_subplot(223)
	#ax3.scatter(f1R550data,f2R550data,  c="green",s=10,marker=".")
	ax3.set_title("R550")
	ax3.set_xlim(0,150)
	ax3.set_ylim(0,150)
	ax3.set_ylabel(sensor2)
	ax3.set_xlabel(sensor1)
	ax3.annotate("r = %.3f" % cor_r550,xy=(0.25, 0.94),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
	
	#regression line
	x_r550 = numpy.array(f1R550data).flatten()
	y_r550 = numpy.array(f2R550data).flatten()
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_r550,y_r550)
	line_r550 = slope*x_r550+intercept
	ax3.plot(x_r550,y_r550,'go', x_r550, line_r550,"b-",linewidth=0.5,markersize=1)
	
	#R1600
	cor_r1600, bias_r1600,RMSD_r1600 = calculateStatistics(f1R1600data,f2R1600data,stat_file,"R1660")
	ax4 = fig.add_subplot(224)
	ax4.set_xlim(0,30)
	ax4.set_ylim(0,30)
	#ax4.scatter(f1R1600data,f2R1600data,  c="red",s=10,marker=".")
	ax4.set_title("R1600")
	ax4.set_xlabel(sensor1)
	ax4.annotate("r = %.3f" % cor_r1600,xy=(0.25, 0.94),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
	
	#regression line
	x_r1600 = numpy.array(f1R1600data).flatten()
	y_r1600 = numpy.array(f2R1600data).flatten()
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_r1600,y_r1600)
	line_r1600 = slope*x_r1600+intercept
	ax4.plot(x_r1600,y_r1600,'ro', x_r1600, line_r1600,"b-",linewidth=0.5,markersize=1)

	
	fig.savefig("%s_ndsi.png" % (out_diff[:-4]))
	
