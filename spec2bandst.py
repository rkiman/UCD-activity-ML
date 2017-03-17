from astropy.io import fits
import math as math
import numpy as np
import glob
from aux_functions import *

def import_indices(file_location):
	'''
	This functions reads the file with the name of the bands and the indices, and returns two lists: one with names of bands and other with the indices. 
	'''
	indices_text = open (file_location, 'r') 

	indices = []
	band_name = []
	
	for line in indices_text:
		if not line.startswith("#"):
			line_text = line.split(', ')
			line_num = []
			for text in line_text[2:]:
				line_num.append(float(text))
			indices.append(line_num)
			band_name.append(line_text[1])
	return band_name,indices
	
def import_fits_data(name):
	'''
	This function opens the fit file. It can be use to analize what is inside the file.
	'''
	hdulist = fits.open(name)
	
	#to get information of the file
	#hdulist.info() #gives number and info of hdu objects in file
	#print(hdulist[1].columns.info)
	#hdulist[1].header['comment'] #gives comment on the hdu selected
	#print(hdulist[1].columns) #if hdulist is data, gives label of columns
	#print(hdulist[1].data[0]) #to access lines
	#print(hdulist[1].data.field(11)) #to access columns 
	#print(hdulist[1].data.shape) #shape of data in hdulist selected
	#print(hdulist[1].data['energy']) #to access columns but with the name 

	return hdulist
	
def get_spectrum(hdulist):
	'''
	This funtion depends on the type of data in the fits file. Analyze fits file and adapt code, if necessary.
	'''
	scidata = hdulist[1].data
	flux = scidata.field(0)
	sigflux = scidata.field(2)
	wavelength = scidata.field(1)
	wavelength_scaled = [10**x for x in wavelength]
	
	return wavelength_scaled,flux,sigflux

def get_avgflux(wavelength_big,flux_big,sigflux_big,startw,endw):
	'''
	Determine the avg flux in the region
	'''	
	
	wavelength_big_array = np.asarray(wavelength_big)   #needs to be an array for next line
	index_range = [i for i,v in enumerate(wavelength_big_array >= startw) if v]
	wavelength = [wavelength_big[x] for x in index_range]
	#I think that calculate the pix scale at the begining of the range is important because the pix_scale change during the spectrum, but not so much in a little range as the one selected here, but it has to be a small range.
	pix_scale = wavelength[1] - wavelength[0]
	
	index_range = [i for i,v in enumerate(np.logical_and(wavelength_big_array+pix_scale/2 >= startw, wavelength_big_array-pix_scale/2 <= endw)) if v]
	
	wavelength = [wavelength_big[x] for x in index_range]
	flux = [flux_big[x] for x in index_range]
	if len(sigflux_big) >= 1:
		sigflux = [sigflux_big[x] for x in index_range]
	
	num_pixels = len(wavelength)

	if num_pixels > 1:
		#determine the fractional pixel value for the pixels on the edge of the region
		frac_ini = (wavelength[0] + pix_scale/2 - startw) / pix_scale
		frac_end = (endw - (wavelength[num_pixels-1] - pix_scale/2)) / pix_scale		
		
		sumflux = 0
		sumsigma = 0

		for i in range(0, num_pixels):
			if i == 0:
				pixflux = frac_ini*flux[i]
			elif i == num_pixels - 1:
				pixflux = frac_end*flux[i]
			else:
				pixflux = flux[i]
			sumflux = sumflux + pixflux

			
			#if len(sigflux) > 1:
			#	if i == 0:
			#		pixsigflux = frac_ini**2*sigflux[i]**2
			#	elif i == num_pixels - 1:
			#		pixsigflux = frac_end**2*sigflux[i]**2
			#	else:
			#		pixsigflux = sigflux[i]**2
			#sumsigma = sumsigma + pixsigflux
			
		
		realpix = num_pixels - 2 + frac_ini + frac_end #fractional pixel value
		avgflux = float(sumflux)/realpix
		
		
		#use the sample variance if the sigma
		#spectrum is not present to estimate uncertainty
		#if len(sigflux) > 1:
		#	avgsigflux = math.sqrt(sumsigma)/realpix 
		#else:
		#	avgsigflux = math.sqrt(sum((flux - mean(flux))^2)/(num_pixels-1))/math.sqrt(num_pixels)	
	
		sigma_aux = [(x - np.mean(flux))**2 for x in flux]
		avgsigflux = math.sqrt(sum(sigma_aux)/(num_pixels)) #put num_pixels instead of num_pixels - 1
	if num_pixels == 1:
		frac = (endw - startw)/pix_scale
		avgflux = frac*flux[0]
		
	return avgflux,avgsigflux

def calculate_bandst_file():
	'''
	This functions creats a file with the bandstreght calculated and the spectral type for each spectrum in the fits files. Just using this function in the main file is enough.
	It is prepared to get the spectral type from a different file.
	'''
	
	band_name,indices = import_indices('/home/ro/Documents/Kelle_Viviana/DML/indices_optical.lis')

	f1 = open('/home/ro/Documents/Kelle_Viviana/bandst.txt','w') 
	f1.write('#' + '\t' + 'file' + '\t' + 'ST' + '\t')

	for name in band_name:
		f1.write(name + '\t' + 'sig' + name + '\t')
	f1.write('\n') 

	files2open = glob.glob('/home/ro/Documents/Kelle_Viviana/SDSS_M_Spectra/*.fits')

	spectraltype = np.loadtxt('/home/ro/Documents/Kelle_Viviana/spectraltype.txt')
	num_bands = len(band_name)
	num_fitsfiles = len(files2open)
				  
	for name in files2open:
		hdulist = import_fits_data(name)
		wavelength, flux, sigflux = get_spectrum(hdulist)

		num = int(filter(str.isdigit, name))
		f1.write(str(num) + '\t')
	
		for x in range(1,len(spectraltype[:,0])):
			if(spectraltype[x,0] == num):
				st = spectraltype[x,1]
		f1.write(str(st) + '\t')	

		for i in range(0,num_bands):
			Fw,sigfluxw = get_avgflux(wavelength,flux,sigflux,indices[i][0],indices[i][1])
			Fs,sigfluxs = get_avgflux(wavelength,flux,sigflux,indices[i][2],indices[i][3])
			bandst = float(Fw)/Fs
			bandst_sig = error(Fw,Fs,sigfluxw,sigfluxs)
			f1.write(str(bandst) + '\t' + str(bandst_sig) + '\t')
				
		f1.write('\n')	

	f1.close() 	
	
	
def create_bandst_file(file_name):
	'''
	This functions creats a file with the bandstreght calculated and the spectral type for each spectrum in the fits files. Just using this function in the main file is enough.
	It is prepared to get the spectral type from a different file.
	'''
	
	f1 = open('/home/ro/Documents/Kelle_Viviana/bandst.txt','w') 
	f1.write('#' + '\t' + 'file' + '\t' + 'ST' + '\t')
	f1.write('TiO5\tsig\tCAH1\tsig\tHalpha\tsig')
	f1.write('\n') 

	hdulist = import_fits_data(file_name)
	scidata = hdulist[1].data
	SPT = scidata.field(11)
	TiO5 = scidata.field(32)
	sigTiO5 = scidata.field(43)
	CAH1 = scidata.field(33)
	sigCAH1 = scidata.field(34)
	Halpha = scidata.field(12)
	sigHalpha = scidata.field(13)
	
	c = zip(SPT,TiO5,sigTiO5,CAH1,sigCAH1,Halpha,sigHalpha)
	f1.write(str(bandst) + '\t' + str(bandst_sig) + '\t')
				
	f1.write('\n')	

	f1.close() 	
