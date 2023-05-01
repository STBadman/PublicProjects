import numpy as np 
import astropy.units as u 
import matplotlib.pyplot as plt  
from astropy.coordinates import SkyCoord,cartesian_to_spherical   
import sunpy.map 
from sunpy.coordinates import frames 
from sunpy.data.sample import AIA_171_IMAGE    
from sunpy.time import parse_time
from scipy.interpolate import griddata
import sunpy.coordinates
import sunpy.coordinates.wcs_utils
from sunpy.net import Fido, attrs as a
from datetime import datetime,timedelta
import warnings,pickle

### Define function to repeat this work with the mapping:
# Input: sunpy.map.Map instance -> [carrington lons,carrington lats, interpolated map data]
def get_carr_map(sunpymap,nlon=720,nlat=360,interp_mode="linear",
				 resample_size=None,crop=60) :
	if not isinstance(sunpymap,sunpy.map.map_factory.GenericMap) : return "Input must be sunpy map instance (sunpy.map.map_factory.GenericMap)"
	
	# Crop to full disk
	map_crop = sunpymap.submap(SkyCoord([-1100, 1100]*u.arcsec,
										[-1100, 1100]*u.arcsec,
										frame=sunpymap.coordinate_frame))

	# If resample_size arg passed (integer), resample map to (resample_size,resample_size)
	if isinstance(resample_size,int) : map_crop=map_crop.resample(u.Quantity([resample_size,resample_size]*u.pixel))
	
	# Get data as ndarray
	map_dat = map_crop.data

	# Map Central Longitude (stonyhurst)
	M0=map_crop.coordinate_frame.observer.lon.value % 360
	
	# Parse obs time to datetime
	date_obs=map_crop.observer_coordinate.frame.obstime.value
	dt_obs = parse_time(date_obs).datetime

	# Generate array of helioprojective coordinates describing pixel locations
	ref_skycoord = map_crop.reference_coordinate
	ref_pix_x,ref_pix_y = (map_crop.reference_pixel.x,
						   map_crop.reference_pixel.y
						  )
	deltx,delty=map_crop.scale.axis1,map_crop.scale.axis2
	xarr = ((np.arange(map_dat.shape[1])*u.pix-ref_pix_x)*deltx 
			+ ref_skycoord.Tx)
	yarr = ((np.arange(map_dat.shape[0])*u.pix-ref_pix_y)*delty 
			+ ref_skycoord.Ty)
	coords_raw_x,coords_raw_y = np.meshgrid(xarr,yarr)
	helioproj_pixel_coords = SkyCoord(coords_raw_x.flatten(),
									  coords_raw_y.flatten(),
									  frame=map_crop.coordinate_frame)

	# Transform pixel coordinates to heliographic stonyhurst and carrington coordinates. 
	# Then Generate list of longitudes and latitudes corresponding to pixels in original data
	# Throws warningsb because of Nan values
	heeq_pixel_coords=helioproj_pixel_coords.transform_to(frames.HeliographicStonyhurst)
	lat_heeq = heeq_pixel_coords.data.lat
	lon_heeq = heeq_pixel_coords.data.lon
	
	### Interpolate irregularly gridded data to a uniform grid in lat and lon

	# Define interpolated grid
	lat_interp = np.linspace(-np.radians(90),np.radians(90),nlat)
	lon_interp = np.radians(np.linspace(0,360,nlon))

	# We do this in stony hurst coordinates but shift the 0 point to 180 degrees to the whole image
	# is in the positive longitude range. We shift to carrington frame at the end.
	lon=((lon_heeq[~np.isnan(lon_heeq)]+np.radians(180)*u.rad) % (np.radians(360)*u.rad)).value
	lat=lat_heeq[~np.isnan(lon_heeq)].value # Remove nan data points (breaks interpolator)
	points = np.array([lon,lat]).T
	values = map_dat.flatten()[~np.isnan(lon_heeq)]
	out_mesh = np.meshgrid(lon_interp,lat_interp)

	#Perform interpolation
	map_interp = griddata(points,values,tuple(out_mesh),method=interp_mode)

	# Restore flattened data to 2d array with dimensions
	# of interpolated grid.
	map_interp=map_interp.reshape(out_mesh[0].shape)  
	
	# Make NaN all data outside +/-60 deg of observation center
	map_interp[:,]
	
	# Correct to carrington coords
	L0 = sunpy.coordinates.sun.L0(dt_obs).value # Carrington coordinate of sub-earth point
	lon_corrected = ((np.degrees(lon_interp)-180)+L0)%360 #shift and enforce [0:360] range
	sort_inds = np.argsort(lon_corrected) # Get indices of lon values when lons are ordered to ascend
	lon_corrected = np.sort(lon_corrected) # order lons in ascending order
	map_interp = map_interp[:,sort_inds]# order EUV data by ascending longitude
	
	# Throw away data greater than 60 deg from lon center
	M0_carr = (M0+L0)%360
	if (M0_carr-crop >= 0) and (M0_carr + crop <= 360) : 
		naninds = np.where((lon_corrected < M0_carr-crop)
							| (lon_corrected > M0_carr+crop)
						  )[0]
	else :
		naninds = np.where((lon_corrected > ((M0_carr+crop) % 360))
							& (lon_corrected < ((M0_carr-crop)%360))
						  )[0]
	map_interp[:,naninds]=np.nan
	
	return lon_corrected,np.degrees(lat_interp),map_interp

def plot_carr_map(lons,lats,values,sunpymaps=None,pcmax=90,pcmin=0,
				  cmap='sdoaia193',figax=None,sinlat=False):
	#Plot
	if sunpymaps is not None : date_obs=sunpymaps[0].observer_coordinate.frame.obstime.value
	if figax is None : fig,ax=plt.subplots(figsize=(20,10))
	else : fig,ax=figax
	if sinlat : 
		ax.pcolormesh(lons,np.sin(np.radians(lats)),values,vmin=np.nanpercentile(values,pcmin),
					   vmax=np.nanpercentile(values,pcmax),cmap=cmap)
		ax.set_aspect(90)
		ax.set_ylabel("Sine(Carrington Latitude)",fontsize=18)
	else : 
		ax.pcolormesh(lons,lats,values,vmin=np.nanpercentile(values,pcmin),
						vmax=np.nanpercentile(values,pcmax),cmap=cmap)
		ax.set_aspect(1)
		ax.set_ylabel("Carrington Latitude / Deg",fontsize=18)
		ax.set_yticks(np.arange(-90,91,30))
	ax.set_xlabel("Carrington Longitude / Deg",fontsize=18)
	
	
	ax.set_xticks(np.arange(0,361,30))
	ax.grid(which="major")

	if sunpymaps is not None :ax.set_title(f"{date_obs} {[(sunpymap.observatory,sunpymap.instrument) for sunpymap in sunpymaps]}",fontsize=18)
	if figax is None : plt.show()

def time_string(dt) : return f"{dt.year}-{dt.month:02d}-{dt.day:02d}T{dt.hour:02d}:{dt.minute:02d}:{dt.second:02d}"

def download_and_gen_carrmap(dt,resample_size=None) :
	# Saves fits files to ~/sunpy/data
	stereoA = (a.vso.Source('STEREO_A') &
			  a.Instrument('EUVI') &
			  a.Time(time_string(dt), time_string(dt+timedelta(minutes=60)))
			 )
	stereoB = (a.vso.Source('STEREO_B') &
			  a.Instrument('EUVI') &
			  a.Time(time_string(dt), time_string(dt+timedelta(minutes=60)))
			 ) 
	aia = (a.Instrument('AIA') &
		   a.Time(time_string(dt), time_string(dt+timedelta(minutes=60)))
		  )
	wave = a.Wavelength(19 * u.nm, 20 * u.nm)
	res = Fido.search(wave, aia | stereoA | stereoB)
	files=[]
	for jj in range(len(res)) :
		try : files.append(Fido.fetch(res[jj,0]))
		except : ""

	maps = [m for m in sunpy.map.Map(files)]
	carrmaps = [get_carr_map(map_,resample_size=resample_size) for map_ in maps]

	lons,lats=carrmaps[0][0:2]

	# Quick process maps to have same floor value (0) and same median
	quick_procs_list=[]
	for arr in carrmaps[:] :
		arr=arr[2]
		quick_proc = arr-np.nanmin(arr)
		quick_procs_list.append(quick_proc/np.nanmedian(quick_proc))
								   
	combined = np.nanmean(np.array(quick_procs_list), axis=0)
	
	return lons,lats,combined,maps    

def get_M0(map_) : return map_.coordinate_frame.observer.lon.value % 360

def front_end(dt,resample_size=512,save=False) :

	lon,lat,combined,maps=download_and_gen_carrmap(dt,resample_size=resample_size)

	fig=plt.figure(figsize=(16,12))
	# Drawing stonyhurst grid on STEREO throws a lot of warnings...
	warnings.filterwarnings('ignore') 

	int_init=201
	int_init += 10*len(maps)
	share=False
	ax1=plt.subplot(int_init)
	for map_ in sorted(maps,key=get_M0):
		ax=plt.subplot(int_init,projection=map_)
		crop = map_#map_.submap(SkyCoord([-1100, 1100]*u.arcsec,
				   #                 [-1100, 1100]*u.arcsec,
				   #                 frame=map_.coordinate_frame))
		pl=crop.plot(cmap='sdoaia193',vmax=2000)
		crop.draw_grid()
		int_init += 1
		
	ax_carr = plt.subplot(212)
	plot_carr_map(lon,lat,combined,maps,pcmax=90,figax=(fig,ax_carr),sinlat=False)
	plt.yticks(np.arange(-90,90,15))
	plt.xticks(np.arange(0,360,15))
	plt.grid()

	plt.tight_layout()

	if not save : plt.show()
	elif isinstance(save,str) : 
		plt.savefig(f"{save}/{str(dt)[0:10]}.png")
		plt.close()
	else: 
		sav_def = "./data/"
		print(f"Saving to default location : {sav_def}") 
		plt.savefig(f"{sav_def}/{str(dt)[0:10]}.png")
		plt.close()

def save_png_over_date_range(dt_start,dt_end,sav_def = "./data/") : 
	warnings.filterwarnings('ignore')
	dt=dt_start
	while dt < dt_end : 
		try : front_end(dt,save=sav_def)
		except : print(f"No map generated for : {dt}") 
		dt +=timedelta(days=1) 

def save_dat_over_date_range(dt_start,dt_end,interval_days=1,sav_def = "./data/") : 
	warnings.filterwarnings('ignore')
	resample_size=512
	dt=dt_start
	datas={}
	while dt < dt_end : 
		try :
			lon,lat,combined,_ = download_and_gen_carrmap(dt,resample_size=resample_size)
			datas[dt]=combined
		except : print(f"No map generated for : {dt}") 
		dt +=timedelta(days=interval_days) 
	filename = f"{sav_def}{str(dt_start)[0:10]}_{str(dt_end)[0:10]}_{interval_days}.pkl"
	pickle.dump((lon,lat,datas),file=open(filename,"wb"))
	return filename

def combine_map(pickle_file,path="./data/",pcmin=5,pcmax=90,figax=None,ret_all=False,
				mode=np.nanmin,plot=False,ret=False,cmap='sdoaia193',sinelat=True) :	
	lon,lat,dat_dict=pickle.load(open(f"{path}{pickle_file}","rb"))
	times = list(dat_dict.keys())
	data_stack=np.array(list(dat_dict.values()))
	combined = mode(data_stack,axis=0)
	if plot :
		if figax is None : figax=plt.subplots(figsize=(20,10))
		plot_carr_map(lon,lat,mode(data_stack,axis=0),pcmin=pcmin,
		              pcmax=pcmax,sinelat=sinelat,
					  cmap=cmap,figax=figax)
		figax[1].set_title(f"{len(times)} maps from {str(times[0])[0:10]} to {str(times[-1])[0:10]}",fontsize=20)
		#plt.show()
	if ret : return lon,lat,combined,list(dat_dict.keys())
	elif ret_all : return lon,lat,data_stack,list(dat_dict.keys())
