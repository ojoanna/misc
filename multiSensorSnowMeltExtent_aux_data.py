#!/usr/bin/env python

'''
Program creates auxiliary datasets (DEM, CLC) for wet snow mapping
'''

# -*- coding: utf-8 -*-

import time, os, optparse, sys
import os
import ogr
import pdb
import copy
import glob
import datetime
import math
import numpy
from pyproj import Proj, transform

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
	print _geom
	return _xmin,_xmax,_ymin,_ymax, _geom

def transform_coordinates(_x1,_y1,_inEPSG,_outEPSG):
	
	_inProj=Proj(init=_inEPSG)
	_outProj=Proj(init=_outEPSG)
	_x2,_y2 = transform(_inProj,_outProj,_x1,_y1)
	
	return _x2,_y2

def select_DEM(poly_LS):

	dem=[]
	if shp_eudem.Contains(poly_LS):
		dem="EUDEM"
	elif srtm_shape.Contains(poly_LS):
		dem="SRTMV4"
	else:
		aster_shape.Contains(poly_LS)
		dem="ASTGTM2"
	print dem

	return dem
		
def create_DEM(dem,corner):
	#corner=poly_LS.GetEnvelope() #xmin,xmax,ymin,ymax
	print corner
	xmin_shp,xmax_shp,ymin_shp,ymax_shp
	xmin=int(corner[0])
	ymin=int(corner[2])
	xmax=int(corner[1])
	ymax=int(corner[3])
	
	print "xmin = ", xmin, "ymin = ", ymin, "xmax = ", xmax, "ymax = ", ymax

	demstring=""
	
	#demtmp="Eudem_N%sE%s_N%sE%s_100m_tmp.tif" % (abs(ymin),abs(xmin),abs(ymax),abs(xmax))
	#pdb.set_trace()
	if os.path.exists(demtmp)==False:
		if dem=="EUDEM":
			nodata=-9999
			DEMDIR="/mnt/e3/rsdb/data/DEM/EU_DEM/EUDEM_Tiles"
			for i in range(ymin,ymax):
				for k in range(xmin,xmax):
					if xmin>=0:
						ori="e"
					else:
						ori="w"
						k=abs(k)
					print i,k
					demfile=glob.glob("%s/n%02d%s%03d.zip" % (DEMDIR,i,ori,k))
					print demfile
					#pdb.set_trace()
					if len(demfile)>0: #os.path.exists(demfile[0])==True:
						cmd="unzip  -q -n %s -d ." % (demfile[0])
						print cmd
						os.system(cmd)
					
						demstring=demstring+" "+demfile[0].split("/")[-1][:-4]+".tif"

		if dem=="SRTMV4":
			DEMDIR="/mnt/e3/rsdb/data/DEM/SRTM_V4"
			xmin2_pos,xmax2_pos,ymin2_pos,ymax2_pos=SRTMv4_select(corner)
			nodata=-32768
			for i in range(xmin2_pos,xmax2_pos+1):
				for k in range(ymax2_pos,ymin2_pos+1):
					demfile=glob.glob("%s/srtm_%02d_%02d.zip" % (DEMDIR,i,k))
					#print demfile
					if os.path.exists(demfile[0])==True:
						cmd="unzip -q -n -j '%s' '%s.tif' -d ." % (demfile[0],demfile[0].split("/")[-1][:-4])
						print cmd
						os.system(cmd)
					demstring=demstring+" "+demfile[0].split("/")[-1][:-4]+".tif"

		
		if dem=="ASTGTM2":
			#DEMDIR="/mnt/e3/rsdb/data/DEM/ASTGTM2"
			DEMDIR="/mnt/e3/rsdb/data/DEM.COPIED/ASTGTM2"
			nodata=-9999
			for i in range(ymin,ymax):
				for k in range(xmin,xmax):
					if xmin>=0:
						ori="E"
					else:
						ori="W"
					
						k=abs(k)
						
					print i,k
					try:
						demfile=glob.glob("%s/ASTGTM2_N%02d%s%03d.zip" % (DEMDIR,i,ori,k))
						print demfile
						cmd="unzip -q -n -j %s '%s/%s_dem.tif' -d ." % (demfile[0],demfile[0].split("/")[-1][:-4],demfile[0].split("/")[-1][:-4])
						#cmd="unzip -q -n -j %s '%s_dem.tif' -d ." % (demfile[0],demfile[0].split("/")[-1][:-4])
						os.system(cmd)
						demstring=demstring+" "+demfile[0].split("/")[-1][:-4]+"_dem.tif"
					except:

						pass
					
		print demstring 
		cmd="gdal_merge.py -n %s -a_nodata -9999 -o %s %s" % (nodata,demtmp,demstring)
		print cmd
		os.system(cmd)
		
		#remove temporary files 
		
		if dem=="EUDEM":
			cmd = "rm n*"
		elif dem=="SRTMV4":
			cmd="rm srtm*"
		elif dem=="ASTGTM2":
			cmd="rm ASTG*"
		cmdexec(cmd)
		
	
		
def SRTMv4_select(corner):
	xmin=math.trunc(corner[0])
	ymin=math.trunc(corner[2])
	xmax=int(math.ceil(corner[1]))
	ymax=int(math.ceil(corner[3]))

	if (ymin-int(numpy.ceil(ymin/10)*10))<5:
		ymin2=int(numpy.ceil(ymin/10)*10)
	else:
		ymin2=int(numpy.ceil(ymin/10)*10)+5
	if (ymax-int(numpy.ceil(ymax/10)*10))>5:
		ymax2=int(numpy.ceil(ymax/10)*10)+5
	else:
		ymax2=int(numpy.ceil(ymax/10)*10)
	if (xmin-int(numpy.ceil(xmin/10)*10))<=5:
		xmin2=int(numpy.ceil(xmin/10)*10)
	else:
		xmin2=int(numpy.ceil(xmin/10)*10)+5
	
	if (xmax-int(numpy.ceil(xmax/10)*10))>=5:
		xmax2=int(numpy.ceil(xmax/10)*10)+5
	else:
		xmax2=int(numpy.ceil(xmax/10)*10)

	ymin2_pos=range(55,-1,-5).index(ymin2)+1
	if ymax==60:
		ymax2_pos=1
	else:
		ymax2_pos=range(55,-1,-5).index(ymax2)+1
	xmin2_pos=range(-180,180,5).index(xmin2)+1
	xmax2_pos=range(-180,180,5).index(xmax2)+1

	#print "xmin=%d xmin2=%d xmin_pos=%d | ymin=%d ymin2=%d ymin_pos=%d | srtm_%02d_%02d.zip" % (xmin,xmin2,xmin2_pos,ymin,ymin2,ymin2_pos,xmin2_pos,ymin2_pos)
	#print "xmax=%d xmax2=%d xmax_pos=%d | ymax=%d ymax2=%d ymax_pos=%d | srtm_%02d_%02d.zip" % (xmax,xmax2,xmax2_pos,ymax,ymax2,ymax2_pos,xmax2_pos,ymax2_pos)

	return xmin2_pos,xmax2_pos,ymin2_pos,ymax2_pos


def cut_DEM(_demtmp,_utm,ul_lon_utm,ul_lat_utm,lr_lon_utm,lr_lat_utm):
	
	#print _demtmp
	##print _utm
	#print ul_lon_utm,ul_lat_utm,lr_lon_utm,lr_lat_utm
	_dem_out="%s.tif" % (_demtmp[:-8])
	demepsg=4326
	
	if os.path.exists(_dem_out)==False:
		cmd="gdalwarp -te %f %f %f %f  -tr 100 100 -s_srs EPSG:%s -t_srs EPSG:326%02d -r bilinear -srcnodata -9999 -dstnodata -9999 %s %s" % (float(ul_lon_utm)-15,float(lr_lat_utm)-15,float(lr_lon_utm)+15,float(ul_lat_utm)+15,demepsg,int(_utm),_demtmp,_dem_out)
		print cmd
		os.system(cmd)

	dem_mdi=_dem_out[:-4]+".mdi"
	if os.path.exists(dem_mdi)==False:
		cmd="dimport -i %s -n dem --itype dem --omittval -9999" % (_dem_out)
		print cmd
		cmdexec(cmd)
		
	cmd = "egeoidm96 -o %s.mdi,a -n egeoid96 --sampling=10,10" % (_dem_out[:-4])
	cmdexec(cmd)
	
	cmd = "ari -m %s.mdi,_dem -a %s.mdi,_egeoid96 -o %s.mdi,a -x M+A -n dem_ellips" % (_dem_out[:-4],_dem_out[:-4],_dem_out[:-4])
	cmdexec(cmd)
	
	cmd="rm %s" % _demtmp
	#cmdexec(cmd)
	
	
def prepareSCM(_track,_auxDir,_clcRaw,_utm,ul_lon_utm,ul_lat_utm,lr_lon_utm,lr_lat_utm):
	
	print "------------------------------------------------------PREPARING LAND COVER MASK---------------------------------------------------------------------------------------------"
	scm_out="SCM_track%s_utm%s_100m_corine.tif" % (_track,_utm)
	if os.path.exists(scm_out)==False:
		cmd = "gdalwarp -s_srs EPSG:4326 -t_srs EPSG:326%s -tr 100 100 -te %s %s %s %s %s %s" % (_utm, float(ul_lon_utm)-15,float(lr_lat_utm)-15,float(lr_lon_utm)+15,float(ul_lat_utm)+15,_clcRaw,scm_out)
		cmdexec(cmd)
	
		#import to .mdi 
		cmd="dimport -i %s" % scm_out
		cmdexec(cmd)
	
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":

	#  populating command line parser with options
	parser = optparse.OptionParser()
	
	parser.add_option("--track",
				action="store", type="string", dest="track", default = None,
				help="Specify the track")
	parser.add_option("--out_dir_sar",
			action="store", type="string", dest="out_dir_sar", default = None,
			help="Specify output directory.")
	parser.add_option("--utm",
			action="store", type="string", dest="utm", default = None,
			help="Specify output directory.")
	
	#config parser
	options, args = parser.parse_args(sys.argv[1:])
	track = options.track
	outdir_sar = options.out_dir_sar
	utm=options.utm
	
	#define directories
	trackDir = "%s/track%s" % (outdir_sar,track)
	auxDir = "%s/track%s/aux_s" % (outdir_sar, track)
	kmlDir = "/mnt/h5/projects/S14SCI_snow/NtoS_track168/kml"
	CLCRAW = "/mnt/h5/projects/S14SCI_snow/NtoS_track168/Europe_N2S/AUX/GlobCover2009v2.3_Corine2012v18.5_Transv14ac_llwgs84_W011N35_E050N72_0.005_cryoland_v3.tif"
	outbound_shp = "%s/AOI_Europe_track%s.shp" % (kmlDir,track)
	print outbound_shp
	
	#calculate boundin box in latlon and utm
	xmin_shp,xmax_shp,ymin_shp,ymax_shp, geom=get_coordinates(outbound_shp)
	print xmin_shp, xmax_shp, ymin_shp,ymax_shp
	inEPSG="epsg:4326"
	outEPSG="epsg:326%s" % utm
	xmin_utm,ymin_utm=transform_coordinates(xmin_shp,ymin_shp,inEPSG,outEPSG)
	
	xmax_utm,ymax_utm=transform_coordinates(xmax_shp,ymax_shp,inEPSG,outEPSG)
	
	print xmin_utm, xmax_utm, ymin_utm,ymax_utm
	
	#DEMs shapefiles
	#srtm
	
	poly_srtm_wkt=("POLYGON((-180. 60.,180. 60.,180. 0.,-180. 0.,-180. 60.))")
	srtm_shape = ogr.CreateGeometryFromWkt(poly_srtm_wkt)
	
	#astgm
	#aster_shape=shapely.geometry.box(-40, 0, 180, 90)
	poly_aster_wkt=("POLYGON((-40. 90.,180. 90.,180. 0.,-40. 0.,-40. 90.))")
	aster_shape = ogr.CreateGeometryFromWkt(poly_aster_wkt)
	
	#EUDEM extent
	eudem_shape="/mnt/h5/projects/S14SCI_snow/NtoS_track168/script/aux/EUDEM_outlines_tiles.shp"
	driver = ogr.GetDriverByName('ESRI Shapefile')
	eudem_dataSource = driver.Open(eudem_shape, 0)
	eudem_layer = eudem_dataSource.GetLayer(0)
	eudem_feature = eudem_layer.GetFeature(0)
	shp_eudem = eudem_feature.GetGeometryRef()
	
	if os.path.exists(trackDir)!=True:
		os.mkdir(trackDir)
	else:
		pass
		
	if os.path.exists(auxDir)!=True:
		os.mkdir(auxDir)
	
	os.chdir(auxDir)
	
	#dem=select_DEM(geom)
	dem="ASTGTM2"
	demtmp="Eudem_N%sE%s_N%sE%s_100m_tmp.tif" % (int(abs(ymin_shp)),int(abs(xmin_shp)),int(abs(ymax_shp)),int(abs(xmax_shp)))
	corner=[xmin_shp,xmax_shp,ymin_shp,ymax_shp]
	create_DEM(dem,corner)

	dem=glob.glob("%s/Eudem*.tif"% auxDir) 
	demtmp=dem[0]
	cut_DEM(demtmp,utm,xmin_utm,ymax_utm,xmax_utm,ymin_utm)
	prepareSCM(track,auxDir,CLCRAW,utm,xmin_utm,ymax_utm,xmax_utm,ymin_utm)
