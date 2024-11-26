import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
import matplotlib.ticker as mticker
from scipy.optimize import fsolve


fontsize=15
tick_dir='out'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size



beta = 2.75
def sigmaHI_over_sigmaHI0_numeric(nu):
    return (nu/nu_HI)**(-beta)
def sigma_over_sigmaHI0_analytic(alpha):
    num = alpha * (1 - 4**(-alpha-beta))
    den = (alpha + beta)*(1 - 4**(-alpha))
    return num/den
def find_root_sigma(alpha,*C):
    return sigma_over_sigmaHI0_analytic(alpha)-C

def get_spectral_index(I_nu_r,sigma_nu,lognu):
    avg_sigma = np.trapz(sigma_nu*I_nu_r,x=lognu)/np.trapz(I_nu_r,x=lognu) # find photon-averaged sigma
    alpha_calc = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma)[0]
    return alpha_calc

def get_tau_avg(r_cMpch,lambda_912): #both cMpc/h
    alpha=1.5
    beta=2.75
    beta_N = 1.6
    #lambda_912 = 16.6 #cMpc/h   #/0.67/(1+z) * cm_per_kpc * kpc_per_Mpc 
    kappa_912 = 1/lambda_912 #h/cMpc
    A = alpha/(alpha+beta)
    B = (1-4**(-alpha-beta*(beta_N-1))) / (1-4**(-alpha)) /(beta_N-1) 
    kappa_avg = kappa_912 * A * B 
    lambda_avg = 1/kappa_avg #cMpc/h
    return kappa_avg*r_cMpch #h/cMpc * cMpc/h = unitless. converting to proper is not necessary

def get_mean_bins(bins,num_bins,I_prime):
    means = []
    # Calculate the mean of I_prime for each bin
    for i in range(1, num_bins + 1):
        # Select the I_prime values that correspond to the current bin
        bin_mask = bins == i
        I_prime_in_bin = I_prime[bin_mask]

        # Calculate the mean of I_prime within the bin
        if len(I_prime_in_bin) > 0:
            mean_value = np.median(I_prime_in_bin)
        else:
            mean_value = np.nan  # Handle empty bins by setting mean as NaN
        # Store the mean value
        means.append(mean_value)
    return means

def power(x, b, c): #power law
    return b * x ** -c

file_path = "./ndot_catalog_062"
X = np.fromfile(file_path,dtype=np.float64)
N_sources = X.shape[0]//5 #5 columns, N_sources rows

# first 3 columns are for x,y,z coordinates. the 4th one is mass. the 5th one is for luminosity 
ndot_catalog = X.reshape((N_sources,5))
# ndot_subcat = ndot_catalog[:1000] #only subset because hard to visualize 2e+5 points at once
Ndot = ndot_catalog[:,4]
log_Ndot = np.log10(Ndot)
a = (log_Ndot-min(log_Ndot)+1e-2)
stdized_log_Ndot = a/max(a) #standardized logNdot

L_box = 40 #cMpc/h
x0,y0,z0 = 12,0.8,32
xc,yc,zc = 20,20,20
dx,dy,dz= xc-x0,yc-y0, zc-z0
#### shifting each coordinate so neutral island is in the center. If <0 or >L, the mod operator takes care of things
x = np.mod(ndot_catalog[:,0]+dx,L_box)
y = np.mod(ndot_catalog[:,1]+dy,L_box)
z = np.mod(ndot_catalog[:,2]+dz,L_box)

distances_cMpch = np.sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2) #now Neutral island is at xc,yc,zc. Simple to find the distance.
distances_pMpc = distances_cMpch /0.68/(1+5.7)

num_bins = 100
min_dist = np.min(distances_cMpch)
max_dist = np.max(distances_cMpch)
binsize = (max_dist - min_dist)/num_bins
min_dist = 1
bin_edges = np.linspace(min_dist, max_dist, num_bins + 1)
bin_centers = np.linspace(min_dist+binsize/2,max_dist-binsize/2,num_bins)

bins = np.digitize(distances_cMpch, bin_edges)
# Initialize a list to store the mean of I_prime within each bin

tau_avg1 = 0
tau_avg2 = get_tau_avg(distances_cMpch,16.6) #kappa_avg*distances_cMpch
tau_avg3 = get_tau_avg(distances_cMpch,16.6-7.5) 
I_prime1 = ndot_catalog[:,4]/(distances_pMpc)**2 * np.exp(-tau_avg1)
I_prime2 = ndot_catalog[:,4]/(distances_pMpc)**2 * np.exp(-tau_avg2)
I_prime3 = ndot_catalog[:,4]/(distances_pMpc)**2 * np.exp(-tau_avg3)

means1 = np.array(get_mean_bins(bins,num_bins,I_prime1))
means2 = np.array(get_mean_bins(bins,num_bins,I_prime2))
means3 = np.array(get_mean_bins(bins,num_bins,I_prime3))

interp1 = np.poly1d(np.polyfit(bin_centers,np.log10(means1),5))
interp2 = np.poly1d(np.polyfit(bin_centers,np.log10(means2),5))
interp3 = np.poly1d(np.polyfit(bin_centers,np.log10(means3),5))
# p = np.poly1d(z)

plt.figure(figsize=(5,5))
# plt.scatter(distances_cMpch,I_prime,s=1,alpha=0.5)
plt.scatter(bin_centers,means1,s=5,label="lambda_912 = infinity",color='red')#"{:d} bins".format(num_bins))
plt.scatter(bin_centers,means2,s=5,label="lambda_912 = 16.6 cMpc/h",color='blue')#"{:d} bins".format(num_bins))
plt.scatter(bin_centers,means3,s=5,label="lambda_912 = 9.1 cMpc/h",color='green')#"{:d} bins".format(num_bins))
plt.plot(bin_centers,10**interp1(bin_centers),color='black',ls="--")
plt.plot(bin_centers,10**interp2(bin_centers),color='black',ls="--")
plt.plot(bin_centers,10**interp3(bin_centers),color='black',ls="--")

plt.xlabel("distance cMpc/h")
plt.ylabel(" Ndot exp[-tau] / r^2 ")
plt.yscale("log")
plt.xlim(0,20)
# plt.axvline(20,color='red')
plt.legend()
# plt.show()
plt.savefig("source_distribution_v2.pdf")
plt.close()

####################################################################################
### Use the interp2 for the weighting function since it has the shortest MFP ###
####################################################################################
Inu_avg = np.loadtxt("Inu_avg.txt")
alpha_r = np.loadtxt("alpha_stack.txt")[0]

N_nu = 30
h_eV = 4.135667e-15
nu_HI = 13.6/h_eV
nu = np.logspace(np.log10(nu_HI),np.log10(4*nu_HI),N_nu) #units of Hz
lognu = np.log(nu)

# R = 5.900
R = 5.5
N_r = 22000
R_cMpc_over_h = R*(1+5.7)*0.68
# N_r = 23600

spatial_fractions = np.array([0,0.2,0.4,0.6,0.8])
travel_dists = spatial_fractions*R_cMpc_over_h
spatial_idxs = np.asarray(spatial_fractions*N_r,int)
systematic_offset = np.max(alpha_r)-1.5
alpha_dists = np.array([alpha_r[i]-systematic_offset for i in spatial_idxs])

cell_size = R_cMpc_over_h/N_r
first_cell = int((1.5)/cell_size)
last_cell = int((20)/cell_size)
N_chunk = last_cell-first_cell
positions = np.linspace(1.5,20,N_chunk) #cMpc/h
source_distribution = interp2(positions)


unnormed_weighted_spectrum = np.sum(np.multiply(Inu_avg[first_cell:last_cell].T,source_distribution),axis=1)
weighted_spectrum = unnormed_weighted_spectrum/np.sum(source_distribution)
popt, pcov = curve_fit(power, nu, weighted_spectrum, p0=[1.0,1.0])
print(popt)


# sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
# avg_sigma = np.trapz(sigma_nu*Inu_avg[0],x=lognu)/np.trapz(Inu_avg[0],x=lognu)
# alpha_calc = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma)[0] # find the root
# print(alpha_calc)

sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
avg_sigma = np.trapz(sigma_nu*weighted_spectrum,x=lognu)/np.trapz(weighted_spectrum,x=lognu)
alpha_calc = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma)[0] - 0.01# find the root
print(alpha_calc)
# I_naught=9e-23
# I_test = popt[0]*(nu)**(-alpha_calc)


fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(10,5))
ax[1].scatter(nu/nu_HI,weighted_spectrum,s=10,color='black',
	# label="average intensity weighted by source distribution")
	label=r"avg. intensity weighted by $\frac{\dot N_\mathrm{sources}}{r^2} e^{-\tau}$")

ax[1].plot(nu/nu_HI,power(nu, *popt),label=r"$\alpha_\mathrm{eff}=$"f"{popt[-1]:.2f}",color='black',ls='dashed')
print(np.abs(power(nu, *popt)-weighted_spectrum)/power(nu, *popt) *100)
# plt.plot(nu/nu_HI,Inu_avg[0],color='black',label="alpha=1.5")
# ax[1].plot(nu/nu_HI,I_test,color='purple')

colors=["blue","orange","green","red","purple"]
p0=[1.0,1.0]
for i in range(len(spatial_idxs)):
	idx = spatial_idxs[i]
	alpha_i = alpha_r[idx]-systematic_offset
	# def power2(x, a, b): #power law
	#     return a + b * x ** -alpha_i
	popt, pcov = curve_fit(power, nu, Inu_avg[idx], p0=p0)
	p0 = popt   
	ax[0].scatter(nu/nu_HI,Inu_avg[idx],s=10,color=colors[i],label=f"{travel_dists[i]:.0f} cMpc/h, "
		+r"$\alpha_{eff}=$"+f"{alpha_i:.2f}")
	ax[0].plot(nu/nu_HI,power(nu,*popt),color=colors[i],ls='--')

# idx = spatial_idxs[1]
# print(idx)
# def power2(x, a, b): #power law
# 	    return a + b * x ** (alpha_r[idx]-0.01)
# popt, pcov = curve_fit(power2, nu, Inu_avg[0], p0=[-1.23067516e-23,9.63707915e-08])
# plt.scatter(nu/nu_HI,power2(nu, *popt),label=f"alpha={alpha_r[idx]-0.01:.2f}",color=colors[0])
# print(power2(nu, *popt))
# I_naught = np.max(1.2e-22)
# I_test = I_naught*(nu/nu_HI)**(-alpha_dists[0])
# ax[0].plot(nu/nu_HI,I_test,color="blue",ls='--')

# I_naught = np.max(1.07e-22)
# I_test = I_naught*(nu/nu_HI)**(-alpha_dists[1])
# ax[0].plot(nu/nu_HI,I_test,color="orange",ls='--')

# I_naught = np.max(6.8e-23)
# I_test = I_naught*(nu/nu_HI)**(-alpha_dists[2])
# ax[0].plot(nu/nu_HI,I_test,color="green",ls='--')

# I_naught = np.max(5.1e-23)
# I_test = I_naught*(nu/nu_HI)**(-alpha_dists[3])
# ax[0].plot(nu/nu_HI,I_test,color="red",ls='--')

# I_naught = np.max(3.75e-23)
# I_test = I_naught*(nu/nu_HI)**(-alpha_dists[4])
# ax[0].plot(nu/nu_HI,I_test,color="purple",ls='--')

ax[0].set_xscale("log")
ax[0].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_yscale("log")


# ax[0].tick_params(top=True,right=True,labelbottom=False,which='both',direction='in')
# ax[1].set_xticks(direction='in')

ax[0].set_ylim(1.e-23,1.5e-22)
ax[1].set_ylim(1.e-23,1.5e-22)
xticks= [1,2,3,4]
xlabels = np.asarray(xticks,str)
# yticks = np.array([2,3,4,5,6,7,8,9,10])*1e-23
# ylabels = np.asarray(yticks*1e23,int)
# ax[1].set_yticks(ticks=yticks)
# ax[1].set_yticklabels(ylabels)
ax[0].set_xticks(ticks=xticks)
ax[0].set_xticklabels(xlabels)
ax[0].tick_params(top=True,right=True,which='both',direction='in')
ax[1].set_xticks(ticks=xticks)
ax[1].set_xticklabels(xlabels)
ax[1].tick_params(top=True,right=True,which='both',direction='in')


# ax[1].tick_params(labelbottom=False,which='both')
# ax.set_yticks(ticks=ticks)
# ax.set_yticklabels(labels)

# ax[1].xaxis.set_minor_formatter(mticker.ScalarFormatter())
# ax[1].xaxis.set_major_formatter(mticker.ScalarFormatter())
# ax[1].xaxis.get_major_formatter().set_scientific(False)
#"average intensity weighted by source distribution"
# ax[0].text(x=0.2,y=0.92,s=r"$I_\nu(r)$",bbox=dict(facecolor='white', alpha=0.8),
# 	transform=ax[0].transAxes,fontsize=15)
# ax[1].text(x=0.2,y=0.92,s=r"$\left<I_\nu\right>_r$",bbox=dict(facecolor='white', alpha=0.8),
# 	transform=ax[1].transAxes,fontsize=15)

	# s=r"$\alpha=${}, {}".format(alpha_list[j],reion_hist_list[j]),
 #        bbox=dict(facecolor='white', alpha=0.8),transform=ax[1][j+2].transAxes,color=color_list[j])
 #    im_list.append(im)
# ax[0].set_ylabel(r"$I_\nu(r)$")
# ax[1].set_ylabel(r"$\left<I_\nu\right>_r$")
# ax[0].set_ylabel(r"intensity $\times 10^{23}$")
# fig.supylabel(r'intensity [erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ Hz$^{-1}$]')
# ax[1].set_ylabel("intensity")
# ax[0].set_xlabel("energy [Ry]")
ax[0].set_xlabel("energy [Ry]")
ax[1].set_xlabel("energy [Ry]")
ax[0].legend(loc='upper right',title=r"$\bar I_\nu(r)$",title_fontsize=15,fontsize=10)
ax[1].legend(loc='upper right',title=r"$\left<\bar I_\nu\right>_r$",title_fontsize=15,fontsize=12)
plt.tight_layout()
plt.savefig("power_law_tests.png",dpi=200)
# plt.show()


# Inu_avg = np.loadtxt("Inu_avg.txt")
# h_eV = 4.135667e-15
# nu_HI = 13.6/h_eV
# beta=2.75
# R = 5.900 #pMpc
# R_cMpc_over_h = R*(1+5.7)*0.68
# N_r = 23600

# #user-input. (hopefully don't have to change this)
# N_nu = 30
# nu = np.logspace(np.log10(nu_HI),np.log10(4*nu_HI),N_nu) #units of Hz
# lognu = np.log(nu)


### scatter plot to check ###
# fig,ax = plt.subplots(1,figsize=(5,5))
# ax.scatter(x=x[:1000],
#            y=z[:1000],
#            #s=stdized_log_Ndot*50,
#            alpha=stdized_log_Ndot[:1000],
#            c=Ndot[:1000],cmap='inferno',norm=matplotlib.colors.LogNorm())
# ax.set_ylim(0,40)
# ax.set_xlim(0,40)
# ax.scatter(20,20,marker='*')
# plt.show()

### histogram of distances to check ###
# plt.figure(figsize=(6,4))
# plt.hist(distances_cMpch)
# plt.axvline(20,color='red')
# plt.xlabel("distance cMpc/h")
# plt.ylabel("number of sources")
# print("There are {:d} sources total.".format(len(distances_cMpch)))
# plt.show()


