import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.signal import find_peaks
import math

T = np.load("data_file.npy")    # iteration results: output data table
np.savetxt("data_file.csv",T)
print(T)
ft = [0,5,10,15,20,25,30]
# alpha = [-15,-10,-5,0,5,10,15]  # alpha angle indexing
# a_new = []
b_new = np.arange(-15,16,1)
print(b_new)
# beta =[-30,-25,-20,-15,-10,-5,0,5,10,15,20]      # beta angle indexing

zi = np.zeros((31,7))
zj = np.zeros((31,7))
zk = np.zeros((31,7))

zin = np.zeros((31,7))
zjn = np.zeros((31,7))
zkn = np.zeros((31,7))
zt = np.zeros((31,7))

for n in range(len(T)):
        x_val = int(T[n,2])
        y_val = int(T[n,0])
        zi[x_val,y_val] = T[n,7+7]*100 # scope approach distance
        zj[x_val,y_val] = T[n,10+7]*180/(math.pi)# scope approach angle
        zk[x_val,y_val] = T[n,13+7]*180/(math.pi)# scope entry angle

# # Trying 3D plot
# X,Y = np.meshgrid(beta,ft)
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.contour3D(beta,ft,np.transpose(zi))
# ax.set_title('Scope Approach Distance [UNIT?]')
# ax.set_xlabel('Scope Central Position Angle (deg)')
# ax.set_ylabel('Fetal Tilt Angle (deg)')
# plt.show()

# Normalizing scores (clean this up if there's time)
# min value in score
I_min = np.amin(zi)
J_min = np.amin(zj)
K_min = np.amin(zk)

# range of scores
I_range = np.ptp(zi)
J_range = np.ptp(zj)
K_range = np.ptp(zk)

# normalize each param score over 0:1
for nn in range(len(T)):
        x_val = int(T[nn,2])
        y_val = int(T[nn,0])
        zin[x_val,y_val] = (zi[x_val,y_val] - I_min) / I_range
        zjn[x_val,y_val] = (zj[x_val,y_val] - J_min) / J_range
        zkn[x_val,y_val] = (zk[x_val,y_val] - K_min) / K_range

# normalized score from 0 - 1
zt = zin - 2*zjn - zkn
z_min = np.amin(zt)
z_range = np.ptp(zt)
zt = (zt - z_min)/z_range

# Optimization! Combining all 3 characteristics.
# normalize them first
X,Y = np.meshgrid(b_new,ft)
fig, ax = plt.subplots()
levels = np.linspace(0,1,11)
PI = ax.contourf(b_new,ft,np.transpose(zt),levels=levels,cmap='viridis')
# PJ = ax.contourf(b_new,ft,,cmap='viridis_r',alpha=0.3)
# PK = ax.contourf(b_new,ft,,cmap='viridis_r',alpha=0.3)
ax.set_title('Optimal Scope Position for \nMax Approach Dist, Min Approach & Min Entry Angle (1:2:1)',fontsize=10)
ax.set_xlabel('Scope Offset (cm)')
ax.set_ylabel('Fetal Tilt Angle (deg)')

cbar = fig.colorbar(PI)
# ax.FaceAlpha = 0.3;
plt.show()


# X,Y = np.meshgrid(beta,ft)
# figi, axi = plt.subplots()
# PI = axi.contourf(beta,ft,np.transpose(zi),cmap='viridis')
# axi.set_title('Scope Approach Distance (cm)')
# axi.set_xlabel('Scope Central Position Angle (deg)')
# axi.set_ylabel('Fetal Tilt Angle (deg)')
# cbar = figi.colorbar(PI)
# axi.FaceAlpha = 0.3;
# # axi.clabel(PI,inline=True,fontsize=10)
#
# figj, axj = plt.subplots()
# PJ = axj.contourf(beta,ft,np.transpose(zj),cmap='viridis_r')
# axj.set_title('Scope Approach Angle (deg)')     #*** check Arion's calc- may be tangent, not normal.
# axj.set_xlabel('Scope Central Position Angle (deg)')
# axj.set_ylabel('Fetal Tilt Angle (deg)')
# cbar = figj.colorbar(PJ)
# axj.FaceAlpha = 0.3;
# # axj.clabel(PJ,inline=True,fontsize=10)
#
# figk, axk = plt.subplots()
# PK = axk.contourf(beta,ft,np.transpose(zk),cmap='viridis_r')
# axk.set_title('Scope Entry Angle (deg)')        #*** check Arion's calc- may be tangent, not normal.
# axk.set_xlabel('Scope Central Position Angle (deg)')
# axk.set_ylabel('Fetal Tilt Angle (deg)')
# cbar = figk.colorbar(PK)
# axk.FaceAlpha = 0.3;
# # axk.clabel(PK,inline=True,fontsize=10)

# plt.show()