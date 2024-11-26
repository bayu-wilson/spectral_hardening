import numpy as np
import matplotlib.pyplot as plt

fontsize=12
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
#plt.rcParams['axes.prop_cycle'] = cycler(color=['k','r','b'])
#plt.rcParams['axes.prop_cycle'] = cycler(ls=['-'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size


#R = 5.900 #pMpc
#N_r = 23600
R = 5.5 #pMpc
N_r = 22000
r_pMpc = np.linspace(0,R,N_r)
r_cMpc_over_h = r_pMpc*(1+5.7)*0.68

GammaHI_r = np.loadtxt("GammaHI_r.txt")
alpha_stack = np.loadtxt("alpha_stack.txt")
#alpha_r = np.loadtxt("alpha_stack.txt")[0]

fig,ax = plt.subplots(ncols=2,nrows=1,figsize=(8,4))


gamma_Gaikwad = 0.224e-12
gamma_Gaikwad_lower = gamma_Gaikwad - 0.112e-12
gamma_Gaikwad_upper = gamma_Gaikwad + 0.223e-12
ax[0].plot(r_cMpc_over_h,GammaHI_r,color='black')
ax[0].set_ylabel(r"$\left<\Gamma_\mathrm{HI}\right>$")
#ax[0].axhline(gamma_Gaikwad_lower)
ax[0].axhline(gamma_Gaikwad,color='red',label='Gaikwad+2023')
#ax[0].axhline(gamma_Gaikwad_upper)
ax[0].set_xlabel("Distance traveled through\n"+r"ionized IGM [h$^{-1}$ cMpc]")

systematic_offset = np.max(alpha_stack[0])-1.5
ax[1].plot(r_cMpc_over_h[:-1],alpha_stack[0]-systematic_offset,color='orange',label=r"$\left< \sigma_\mathrm{HI} \right>_\nu = \left< \sigma_\mathrm{HI} \right>_\nu^\mathrm{analytic}(\alpha_\mathrm{eff})$")
#systematic_offset = np.max(alpha_stack[1])-1.5
#ax[1].plot(r_cMpc_over_h[:-1],alpha_stack[1]-systematic_offset,color='green',label="heat weighted")
systematic_offset = np.nanmax(alpha_stack[2])-1.5
ax[1].plot(r_cMpc_over_h[:-1],alpha_stack[2]-systematic_offset,color='blue',label="power law fit")
print(alpha_stack[2])

#ax[1].set_ylim(0.4,1.6)
#ax[1].set_ylim(-0.5,1.6)
#ax[1].set_yticks([0.5, 0.7, 0.9, 1.1, 1.3, 1.5]) 
ax[1].set_ylabel(r"Effective spectral index, $\alpha_\mathrm{eff}$")
ax[1].set_xlabel("Distance traveled through\n"+r"ionized IGM [h$^{-1}$ cMpc]")

for i in range(2):
	ax[i].set_xlim(0,20)#np.max(r_cMpc_over_h))

ax[0].legend()
ax[1].legend()

plt.tight_layout()
plt.savefig("GammaHI_alpha.png",dpi=400)
