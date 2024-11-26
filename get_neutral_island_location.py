import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def recenter_line(X,xx,Ngas):
    X_shift = np.zeros_like(X)
    recenter = int(xx*Ngas)
    for ii in range(Ngas):
        X_shift[ii] = X[ii-recenter]
    return X_shift


file_path = "./ndot_catalog_062"
X = np.fromfile(file_path,dtype=np.float64)
N_sources = X.shape[0]//5 #5 columns, N_sources rows

# first 3 columns are for x,y,z coordinates. the 4th one is mass. the 5th one is for luminosity 
ndot_catalog = X.reshape((N_sources,5)) 
ndot_subcat = ndot_catalog[:1000] #only subset because hard to visualize 2e+5 points at once
Ndot = ndot_subcat[:,4]
log_Ndot = np.log10(Ndot)

stdized_log_Ndot = (log_Ndot-min(log_Ndot)+1e-2) #avoids dividing by zero
stdized_log_Ndot = stdized_log_Ndot/max(stdized_log_Ndot) #standardized logNdot

#############################
### make 3d scatter plot! ###
#############################
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(projection='3d')
sc = ax.scatter3D(ndot_subcat[:,0],ndot_subcat[:,1], ndot_subcat[:,2], 
                  s=stdized_log_Ndot*50,
                  alpha=(stdized_log_Ndot>0.0)*stdized_log_Ndot,
                  c=Ndot,cmap='jet',
                  norm=matplotlib.colors.LogNorm())
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('z-axis')
cbar = fig.colorbar(sc, ax=ax, label='N dot')
#plt.show()
plt.savefig("3D_scatterplot.png")
plt.close()

################################################################################
######### show neutral island  location by finding low GammaHI region ##########
################################################################################
gas_file = "gas_z=05.8349_N96_for_calibration"
Ngas = 96
gas = np.fromfile(gas_file,dtype=np.float32)[1:].reshape((23,Ngas,Ngas,Ngas),order='F')

image = np.sum(gas[8],axis=1)
# image = recenter_image(image,8/40,28/40,Ngas)
# image = recenter_image(image,25/95,80/96,Ngas)
fig,ax = plt.subplots(1,figsize=(5,5))
ax.imshow(image,origin='lower')
ax.scatter(Ngas/2,Ngas/2,color='red',marker='x')

xc,zc = Ngas/2,Ngas/2
x_NI,z_NI = xc-8/40*Ngas, zc+12/40*Ngas
ax.scatter(x_NI,z_NI,marker='*',color='white',label=f"(x,y)=({x_NI/Ngas*40:.1f},{z_NI/Ngas*40:.1f})")
           # Ngas/2+28/40*Ngas %40,marker='*',color='white')
# print(f"The neutral island is located at (x,z)=({x_NI/Ngas*40:.1f},{y_NI/Ngas*40:.1f}) cMph/h")
ticks= np.linspace(0,Ngas,5)
labels = np.asarray(np.linspace(0,40,5),int)
ax.set_xticks(ticks=ticks)
ax.set_xticklabels(labels)
ax.set_yticks(ticks=ticks)
ax.set_yticklabels(labels)
ax.legend()
# plt.show()

########################################
##### overlay the source positions #####
########################################

# fig,ax = plt.subplots(1,figsize=(5,5))
xpoints = (ndot_subcat[:,0])
xpoints = xpoints*Ngas/40
zpoints = (ndot_subcat[:,2])
zpoints = zpoints*Ngas/40
ax.scatter(x=xpoints,
           y=zpoints,
           #s=stdized_log_Ndot*50,
           s=10, color='black',
           alpha=stdized_log_Ndot)
           #c=Ndot#,cmap='inferno',norm=matplotlib.colors.LogNorm())
# ax.scatter(x=xpoints,
#            y=ypoints,
#            #s=stdized_log_Ndot*50,
#            alpha=stdized_log_Ndot,
#            c=Ndot,cmap='inferno',norm=matplotlib.colors.LogNorm())
ax.set_ylim(0,Ngas)
ax.set_xlim(0,Ngas)
plt.savefig("map_GammaHI_with_source_positions_.png")
plt.close()
#plt.show()

gammaHI_box = gas[8]
# x,y,z = int(12/40*Ngas+0.5),int(0.52*Ngas+0.5),int(32/40*Ngas+0.5)
y_NI =int(0.02*Ngas+0.5) #0.52-0.5
# print(f"y={y_cMpc:.1f}",)
print(f"The neutral island is located at (x,y,z)=({x_NI/Ngas*40:.1f},{y_NI/Ngas*40:.1f},{z_NI/Ngas*40:.1f}) cMph/h")
print(f"The neutral island is located at (x,y,z)=({x_NI:.1f},{y_NI:.1f},{z_NI:.1f}) cell")

qwer = recenter_line(gammaHI_box[29,:,77],0.52,Ngas) #x,y,z
plt.plot(qwer)
# plt.plot(gammaHI_box[29,:,77]) #x,y,z)
plt.ylabel("Gamma HI")
plt.axvline(Ngas/2)
plt.savefig("sightline.png")
plt.close()
#plt.show()

### LOCATION OF NEUTRAL ISLAND (x,y,z)=(12.0,0.8,32.0) cMph/h ###
### LOCATION OF NEUTRAL ISLAND (x,y,z)=(12.0,20.8,32.0) cMph/h ###
