import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse_single_table
import matplotlib.image as image
import matplotlib.colors
from mw_plot import MWSkyMap
from mw_plot import MWFaceOn
from astropy import units as u
import astropy.coordinates as apycoords
from matplotlib.lines import Line2D
markers = Line2D.markers

from matplotlib import pyplot as plt
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('Dark2').colors)


plt.rc('font', family='serif')
plt.rc('xtick', labelsize=25)
plt.rc('ytick', labelsize=25)

plot_mw_instance = MWFaceOn(radius=12 * u.kpc, unit=u.kpc, coord='galactic',annotation=True)
plot_mw_instance.title = None  # plot title, or it can be None to show no title
plot_mw_instance.fontsize = 35  # fontsize for matplotlib plotting
plot_mw_instance.figsize = (20, 20)  # figsize for matplotlib plotting
plot_mw_instance.dpi = 100  # dpi for matplotlib plotting
plot_mw_instance.cmap = 'plasma'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
plot_mw_instance.clim = (0,10000) # colorbar range
plot_mw_instance.imalpha = 0.5  # alpha value for the milkyway image
plot_mw_instance.s = 100.0  # make the scatter points bigger
plot_mw_instance.tight_layout = True # whether plt.tight_layout() will be run
plot_mw_instance.initialize_mwplot()

table = parse_single_table("PS_2025.04.09_11.15.50.votable")

methods = table.array['discoverymethod'].data
unique_method = np.unique( table.array['discoverymethod'].data,return_counts=True) [0]

fig,ax = plt.subplots(1,2)

ax[0].axvline(1,lw=30,c='dodgerblue',alpha=0.5)
ax[0].text(0.8,0.5, 'SNOW LINE', color='dodgerblue',fontsize=25,rotation=90)

ax[1].imshow(plot_mw_instance.bg_img, zorder=0, extent=plot_mw_instance.bg_img_ext,  rasterized=True,alpha=0.75)
ratio = []

for index,met in  enumerate(unique_method):

	mask = methods==met

	ax[0].scatter(table.array['pl_orbsmax'].data[mask]/(2.7*table.array['st_mass'].data[mask]),table.array['pl_bmassj'].data[mask]*317.6,s=100,marker = list(markers.keys())[index],alpha=0.75)#,label = met)

	mask = (methods == met) & (~np.isnan(table.array['sy_dist'].data))	
	
	c_dr1 = apycoords.SkyCoord(table.array['ra'].data[mask]*u.deg,table.array['dec'].data[mask]*u.deg, distance=table.array['sy_dist'].data[mask]*u.pc, frame='icrs')

	ax[1].scatter(-c_dr1.galactic.cartesian.x/1000, c_dr1.galactic.cartesian.y/1000,s=100,marker = list(markers.keys())[index],alpha=1,label = met)

	if 'Microlensing' in met:
		ratio.append(table.array['pl_bmassj'].data[mask]*317.6/table.array['st_mass'].data[mask])
		
###OMEGA
#Gaia20bof, Gaia22dkv,AT2021uey, OB231060, OB240034,KB22462

coords = [[184.61816,-63.61816,1.5,np.inf,np.inf],[151.76898,-66.18089,1.2,1.6,1],[324.54500,26.46633,1,4.0,1.3],[270.66679,-35.70483,6,2,1],[266.04333,-39.11322,5,1.8,0.5],[270.79882,-27.48121,4,np.inf,np.inf]]

for coo in coords:
    c_dr1 = apycoords.SkyCoord(coo[0]*u.deg,coo[1]*u.deg, distance=coo[2]*1000*u.pc, frame='icrs')

    ax[1].scatter(-c_dr1.galactic.cartesian.x/1000, c_dr1.galactic.cartesian.y/1000,s=500,marker = '*',c='r',alpha=1)
ax[1].scatter(-c_dr1.galactic.cartesian.x/1000, c_dr1.galactic.cartesian.y/1000,s=500,marker = '*',c='r',alpha=1,label='OMEGA')
    
ax[0].set_xlabel(r'$a [a_{snow}]$',fontsize=35)
ax[0].set_ylabel(r'$M_{planet} [M_\oplus]$',fontsize=35)
ax[0].tick_params(axis='both', pad=10)

ax[0].set_xlim([8*10**-3,10**3])
ax[0].set_ylim([10**-1,10**4])
ax[0].set_yscale('log')
ax[0].set_xscale('log')

ax[1].set_xlabel(r'$X [kpc]$',fontsize=35)
ax[1].set_ylabel(r'$Y [kpc]$',fontsize=35)
ax[1].tick_params(axis='both', pad=10)


ax[1].legend(bbox_to_anchor=(-1.2, 1.0, 2, .102), loc=3, ncol=3, mode="expand", borderaxespad=3,fontsize=15,markerscale=1,frameon=False)


plt.show()


