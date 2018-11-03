#!/usr/bin/python

import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
import astropy
from astropy.cosmology import Planck15
import matplotlib.mlab as mlab 
from matplotlib.colors import LogNorm



'''Plot the flux energy + density phase space for 
	boxfit lightcurves with flux density contours overlaid.

	Currently tailored quite heavily for my use-case...

	Assumes the following file name format:
	lightcurve_6ghz_1e50erg_10cm_30deg.txt

	Currently explicitly assumes model files are for
	z = 0.3, but corrects output flux to source redshift.
	Can change this by editing 'luminosity' 
	variable accordingly. 

	Usage: boxfit_phase_space.py -z redshift -t -time
	
	-t	--time		Specify source time after explosion in days.
	-z	--redshift	Specify source redshift	

	'''



# Read in and process command line options

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Specify directory containing boxfit .txt files', dest='data', type=str, default='None',action='store',required=True)
parser.add_argument('-t', '--time', help='Observer time after explosion in days', dest='time', type=float, default='None', action='store', required=True)
parser.add_argument('-z', '--redshift', help='Source redshift', dest='redshift', type=float, default='None', action='store', required=True)

args = parser.parse_args()
data_direc = args.data
source_time = args.time
source_redshift = args.redshift

# Convert source time to seconds
value = float(source_time * 86400.)


final_result = []
for filename in os.listdir(data_direc):
    if not filename.startswith('.'):

        # Grab energy and density for a given model
        file_energy = filename.split('_')[2]
        file_density = filename.split('_')[3]
        
    
        with open(data_direc+filename) as f:
            
                
            data = f.read() # Read entire file
            data_lines = data.splitlines() # Go line by line
            short_list = data_lines[97:] # Skip header 
            
            rows = []
            for line in short_list:
                rows.append(float(line.split(',')[1])) # Split row into columns
    
            # Find time that is closest to source time
            idx = min(range(len(rows)), key=lambda i: abs(rows[i]-value))
            source_time = short_list[idx].split(',')[1] 
            source_flux = short_list[idx].split(',')[3] # boxfit model predicted flux
            
            search_lines = f.readlines()
        
        with open(data_direc+filename) as f:
            for line in f:
            	# Grab model luminosity distance from boxfit file
                if "luminosity distance" in line:
                    dl = float(line.split(":")[1].split(' ')[1])
            
            # Calculate luminosity of source for model flux (z=0.3)
            luminosity = (1./(1.+0.3)) * 4.*np.pi * dl**2. * float(source_flux) * 1e-3 * 1e-23 
            
            # Calculate flux (in mJy) of source correcting to actual source redshift
            new_flux = (1e23 * 1e3 * luminosity)/(4.*np.pi*Planck15.luminosity_distance(source_redshift).cgs.value**2.* (1.+source_redshift))
                
        # Create row for output array
        # Model energy, model density, source time after explosion, model flux at source time, luminosity distance, flux at new redshift     
        new_row = [file_energy.strip('erg'),file_density.strip('cm'),source_time,source_flux,dl,luminosity,new_flux]
            #print(new_row)
        final_result.append(new_row)

        #final_result = np.array(final_result)
        


#densities_array = [1e-3,1e-2,1e-1,1,10,100]

def plot_phase_space(final_result):


	final_data = np.array(final_result)
	get_energies = final_data[:,0].astype(np.float) # Grab model energies 
	get_densities = final_data[:,1].astype(np.float) # Grab model densities
	get_new_flux = final_data[:,6].astype(np.float) # Grab flux for source redshift


	# Need to grid irregularly spaced data
	ny, nx = 700,700 # Grid size

	# Generate a regular grid to interpolate the data.
	xi = np.linspace(np.min(get_energies), np.max(get_energies), nx)
	yi = np.linspace(np.min(get_densities), np.max(get_densities), ny)
	xi, yi = np.meshgrid(xi, yi)


	# Normalize data
	def normalize_x(data):
		data = data.astype(np.float)
		return (data - np.min(data)) / (np.max(data) - np.min(data))

	def normalize_y(data):
		data = data.astype(np.float)
		return (data - np.min(data)) / (np.max(data) - np.min(data))

	x_new, xi_new = normalize_x(get_energies), normalize_x(xi)
	y_new, yi_new = normalize_y(get_densities), normalize_y(yi)	
	# Linear interpolation (think Delaunay triangulation is preferrable but requires some package...)
	zi = mlab.griddata(x_new, y_new, get_new_flux, xi_new, yi_new,interp='linear')	
	

	# Plot the results
	plt.figure(figsize=(10,8))	
	plt.pcolormesh(xi,yi,zi)
	lvls=[1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100]
	plt.contourf(xi,yi,zi,norm=LogNorm(),levels=lvls)
	
	plt.colorbar(label='Flux [mJy]')
	#plt.clim(1e-3,2)
	
	plt.xlim(1e53,1e55)
	plt.ylim(1e-3,100)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$\rm E_{iso} \ [erg]$')
	plt.ylabel(r'$\rm n \ [cm^{-3}]$')

	plt.tick_params(which='minor',direction='in')
	plt.tick_params(which='major',direction='in')
	plt.tick_params(axis = 'both', which = 'major')
	plt.tick_params(axis = 'both', which = 'minor')
	plt.tick_params(top=True, right=True, which='both') 



	plt.show()

plot_phase_space(final_result)



        
