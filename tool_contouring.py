import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.signal import find_peaks
import math

T = np.load("data_file.npy")    # iteration results: output data table
ft = [0,5,10,15,20,25,30]
# alpha = [-15,-10,-5,0,5,10,15]  # alpha angle indexing
# alpha = [0,5,10,15]
# beta =[-30,-25,-20,-15,-10,-5,0,5,10,15,20]      # beta angle indexing
a_new = np.arange(0,11,1)
b_new = np.arange(-15,16,1)
T[:,0] = 5*T[:,0]

zT = np.zeros((31, 11))

for m in ft:
        zi = np.zeros((31, 11))
        zj = np.zeros((31, 11))
        zk = np.zeros((31, 11))

        zin = np.zeros((31, 11))
        zjn = np.zeros((31, 11))
        zkn = np.zeros((31, 11))
        zt = np.zeros((31, 11))

        for n in range(len(T)):
                x_val = int(T[n,6])
                y_val = int(T[n,5])
                if int(T[n,0]) == m:
                        zi[x_val, y_val-10] = 100*T[n,9+7]  # tool approach distance
                        zj[x_val, y_val-10] = 180*T[n,12+7]/(math.pi)  # tool approach angle
                        zk[x_val, y_val-10] = 180*T[n,15+7]/(math.pi)  # tool entry angle
                        if(zi[x_val,0]) == 0: zi[x_val, 0] = 100 * T[n, 7+7]  # scope approach distance when at midline
                        if(zj[x_val,0]) == 0: zj[x_val, 0] = 180 * T[n, 10+7] / (math.pi)  # scope approach angle when at midline
                        if(zk[x_val,0]) == 0: zk[x_val, 0] = 180 * T[n, 13+7] / (math.pi)  # scope entry angle when at midline

        print(zi)
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
        zin = (zi  - I_min) / I_range
        zjn = (zj - J_min) / J_range
        zkn = (zk - K_min) / K_range

        # normalized score from 0 - 1
        zt = zin - 2*zjn - zkn
        z_min = np.amin(zt)
        z_range = np.ptp(zt)
        zt = np.divide((zt - z_min),z_range)

        # Optimization! Combining all 3 characteristics.
        # normalize them first
        # X, Y = np.meshgrid(beta, ft)
        fig, ax = plt.subplots()
        levels = np.linspace(0, 1, 11)
        PI = ax.contourf(b_new, a_new, np.transpose(zt), levels=levels, cmap='viridis')
        # PJ = ax.contourf(beta,ft,,cmap='viridis_r',alpha=0.3)
        # PK = ax.contourf(beta,ft,,cmap='viridis_r',alpha=0.3)
        ax.set_title('Optimal Tool Position for \nMax Approach Dist, Min Approach & Min Entry Angle (1:2:1); \nFetal Tilt ='+str(m)+'deg',fontsize=10)
        ax.set_xlabel('Scope Offset (cm)')
        ax.set_ylabel('Tool Offset Position (cm)')

        cbar = fig.colorbar(PI)
        plt.savefig(str(m)+'Tool-121-L.png')
        plt.close()

        ## Plot indiv. parameter results
        # print(zi)
        # zT = zi
        # X,Y = np.meshgrid(beta,ft)
        # figi, axi = plt.subplots()
        # PI = axi.contourf(beta,alpha,np.transpose(zi),cmap='viridis')
        # axi.set_title('Tool Approach Distance (cm)'+'     Fetal Tilt = '+str(m)+'deg')
        # axi.set_xlabel('Scope Central Position Angle (deg)')
        # axi.set_ylabel('Tool Offset Position Angle (deg)')
        # cbar = figi.colorbar(PI)
        # plt.savefig(str(m)+'I.png')
        # plt.close()
        #
        # figj, axj = plt.subplots()
        # PJ = axj.contourf(beta,alpha,np.transpose(zj),cmap='viridis_r')
        # axj.set_title('Tool Approach Angle (deg)'+'     Fetal Tilt = '+str(m)+'deg')
        # axj.set_xlabel('Scope Central Position Angle (deg)')
        # axj.set_ylabel('Tool Offset Position Angle (deg)')
        # cbar = figj.colorbar(PJ)
        # plt.savefig(str(m)+'J.png')
        # plt.close()
        #
        # figk, axk = plt.subplots()
        # PK = axk.contourf(beta,alpha,np.transpose(zk),cmap='viridis_r')
        # axk.set_title('Tool Entry Angle (deg)'+'     Fetal Tilt = '+str(m)+'deg')
        # axk.set_xlabel('Tool Central Position Angle (deg)')
        # axk.set_ylabel('Tool Offset Position Angle (deg)')
        # cbar = figk.colorbar(PK)
        # plt.savefig(str(m)+'K.png')
        # plt.close()

# figT = plt.figure()
# axT = plt.axes(projection='3d')
# PT = axT.contour3D(beta,alpha,np.transpose(zT),ft,cmap='viridis_r')
# axT.set_title('Tool Approach Angle (deg)'+'     Fetal Tilt = '+str(m)+'deg')
# axT.set_xlabel('Scope Central Position Angle (deg)')
# axT.set_ylabel('Tool Offset Position Angle (deg)')
# cbar = figT.colorbar(PT)
# plt.show()
# # fig = plt.figure()
# # ax = plt.axes(projection='3d')
# # ax.contour3D(beta,ft,np.transpose(zi))