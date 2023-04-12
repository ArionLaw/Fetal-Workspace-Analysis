from vedo import*
import numpy as np
import pyvista as pv
from math import sin, cos, pi
import os
import sys

CurrentDirectory = sys.path[0]

### For PSM End Effector ###
from psmFK import*
"""
q1 Outer Yaw
q2 Outer Pitch
q3 Insertion
q4 Wrist Roll
q5 Wrist Pitch
q6 Wrist Yaw
"""
# PSM joint max limits: Outer Yaw, Outer Pitch, Insertion, Wrist Roll, Wrist Pitch, Wrist Yaw 
q_max_PSM = [90*pi/180 , 54*pi/180 , 0.240 , 180*pi/180 , 180*pi/180 , 180*pi/180 ] # in radians and meters
q_min_PSM = [-90*pi/180 , -54*pi/180 , 0.000 , -180*pi/180 , -180*pi/180 , -180*pi/180 ] # in radians and meters

# PSM joint soft limits: Outer Yaw, Outer Pitch, Insertion, Wrist Roll, Wrist Pitch, Wrist Yaw 
qlim_u_PSM = [30*pi/180 , 30*pi/180 , 0.2 , 120*pi/180 , 90*pi/180 , 90*pi/180 ] # in radians and meters
qlim_l_PSM = [-30*pi/180 , -30*pi/180 , 0.02 , -120*pi/180 , -90*pi/180 , -90*pi/180 ] # in radians and meters

revolute1_divisions = 5
revolute2_divisions = 5
prismatic_divisions = round((qlim_u_PSM[2] - qlim_l_PSM[2])/0.005)+1
q1 = np.linspace(qlim_l_PSM[0],qlim_u_PSM[0],revolute1_divisions)
q2 = np.linspace(qlim_l_PSM[1],qlim_u_PSM[1],revolute2_divisions)
q3 = np.linspace(qlim_l_PSM[2],qlim_u_PSM[2],prismatic_divisions)

i = 0
j = 0
h = 0
k = 0
PSM_EE_xyz = [[]]
while j in range(prismatic_divisions):
    while k in range(revolute2_divisions):
        while h in range(revolute1_divisions):
            joint_pos = [q1[h] , q2[k] , q3[j] , 0 , 0 , 0]
            base_frame_pos = np.mat([[0],[0],[0],[1]])
            PSM_EE_pos = np.matmul(compute_FK(joint_pos),base_frame_pos)
            PSM_EE_xyz.append([float(PSM_EE_pos[0]),float(PSM_EE_pos[1]),float(PSM_EE_pos[2])])
            i = i + 1

            if j==0 or j==prismatic_divisions-1 or k==0 or k==revolute2_divisions-1 or h==revolute1_divisions-1: 
                h = h + 1
            else:
                h = revolute1_divisions-1
        h = 0
        k = k + 1

    print("Insertion [m]: " , q3[j] )
    k = 0
    j = j + 1
    if abs(PSM_EE_xyz[i][0] - PSM_EE_xyz[i-1][0]) > 0.005:
        revolute1_divisions = revolute1_divisions + 2
        q1 = np.linspace(qlim_l_PSM[0],qlim_u_PSM[0],revolute1_divisions)
        revolute2_divisions = revolute2_divisions + 2
        q2 = np.linspace(qlim_l_PSM[1],qlim_u_PSM[1],revolute2_divisions)

### PSM Linkages ###
revolute1_divisions = 15
revolute2_divisions = 15
prismatic_divisions = round((qlim_u_PSM[2] - qlim_l_PSM[2])/0.005)+1
q1 = np.linspace(qlim_l_PSM[0],qlim_u_PSM[0],revolute1_divisions)
q2 = np.linspace(qlim_l_PSM[1],qlim_u_PSM[1],revolute2_divisions)
q3 = np.linspace(qlim_l_PSM[2],qlim_u_PSM[2],prismatic_divisions)

i = 0
j = prismatic_divisions-1
h = 0
k = 0
PSM_SterileAdapter_xyz = [[]]
while j in range(prismatic_divisions):
    while k in range(revolute2_divisions):
        while h in range(revolute1_divisions):
            joint_pos = [q1[h] , q2[k] , q3[j]]
            base_frame_pos = np.mat([[0],[0],[0],[1]])
            PSM_SterileAdapter_pos = np.matmul(compute_FK(joint_pos),base_frame_pos)
            PSM_SterileAdapter_xyz.append([float(PSM_SterileAdapter_pos[0]),float(PSM_SterileAdapter_pos[1]),float(PSM_SterileAdapter_pos[2])])
            i = i + 1

            if j==0 or j==prismatic_divisions-1 or k==0 or k==revolute2_divisions-1 or h==revolute1_divisions-1: 
                h = h + 1
            else:
                h = revolute1_divisions-1
        h = 0
        k = k + 1

    print("Insertion [m]: " , q3[j] )
    k = 0
    j = j - 1
    if abs(PSM_SterileAdapter_xyz[i][0] - PSM_SterileAdapter_xyz[i-1][0]) > 0.005:
        revolute1_divisions = revolute1_divisions + 2
        q1 = np.linspace(qlim_l_PSM[0],qlim_u_PSM[0],revolute1_divisions)
        revolute2_divisions = revolute2_divisions + 2
        q2 = np.linspace(qlim_l_PSM[1],qlim_u_PSM[1],revolute2_divisions)

### for ECM ###
from ecmFK import*
"""
q1 Outer Yaw
q2 Outer Pitch
q3 Insertion
q4 Roll
"""

# ECM joint max limits: Outer Yaw, Outer Pitch, Insertion, Roll
q_max_ECM = [90*pi/180 , 66*pi/180 , 0.240 , 90*pi/180] # in radians and meters
q_min_ECM = [-90*pi/180 , -44*pi/180 , 0.000 , -90*pi/180] # in radians and meters

# ECM joint soft limits: Outer Yaw, Outer Pitch, Insertion, Roll
qlim_u_ECM = [30*pi/180 , 30*pi/180 , 0.2 , 90*pi/180] # in radians and meters
qlim_l_ECM = [-30*pi/180 , -30*pi/180 , 0.02 , -90*pi/180] # in radians and meters

revolute1_divisions = 15
revolute2_divisions = 15
prismatic_divisions = round((qlim_u_ECM[2] - qlim_l_ECM[2])/0.005)+1
q1 = np.linspace(qlim_l_ECM[0],qlim_u_ECM[0],revolute1_divisions)
q2 = np.linspace(qlim_l_ECM[1],qlim_u_ECM[1],revolute2_divisions)
q3 = np.linspace(qlim_l_ECM[2],qlim_u_ECM[2],prismatic_divisions)

i = 0
j = prismatic_divisions-1
h = 0
k = 0
ECM_xyz = [[]]
while j in range(prismatic_divisions):
    while k in range(revolute2_divisions):
        while h in range(revolute1_divisions):
            joint_pos = [q1[h] , q2[k] , q3[j]]
            base_frame_pos = np.mat([[0],[0],[0],[1]])
            ECM_pos = np.matmul(compute_FK(joint_pos),base_frame_pos)
            ECM_xyz.append([float(ECM_pos[0]),float(ECM_pos[1]),float(ECM_pos[2])])
            i = i + 1

            if j==0 or j==prismatic_divisions-1 or k==0 or k==revolute2_divisions-1 or h==revolute1_divisions-1: 
                h = h + 1
            else:
                h = revolute1_divisions-1
        h = 0
        k = k + 1

    print("Insertion [m]: " , q3[j] )
    k = 0
    j = j - 1
    if abs(ECM_xyz[i][0] - ECM_xyz[i-1][0]) > 0.005:
        revolute1_divisions = revolute1_divisions + 2
        q1 = np.linspace(qlim_l_ECM[0],qlim_u_ECM[0],revolute1_divisions)
        revolute2_divisions = revolute2_divisions + 2
        q2 = np.linspace(qlim_l_ECM[1],qlim_u_ECM[1],revolute2_divisions)      

### Meshing and Volume Visualizations ###
PSM_EE_xyz.pop(0)
PSM_SterileAdapter_xyz.pop(0)
ECM_xyz.pop(0)

mesh_algo_radius = 1
cloudRWS = pv.PolyData(PSM_EE_xyz)
#cloudRWS.plot()
volumeRWS = cloudRWS.delaunay_3d(alpha = mesh_algo_radius)
shellRWS = volumeRWS.extract_geometry()
#shellRWS.plot()

PSM1_RWS = Mesh(shellRWS, c="red", alpha=0.4)
PSM2_RWS = Mesh(shellRWS, c="green", alpha=0.4)

mesh_algo_radius = 1
cloudSterileAdapter = pv.PolyData(PSM_SterileAdapter_xyz)
#cloudSterileAdapter.plot()
volumeSterileAdapter = cloudSterileAdapter.delaunay_3d(alpha = mesh_algo_radius)
shellSterileAdapter = volumeSterileAdapter.extract_geometry()
#shellSterileAdapter.plot()

PSM1_SterileAdapter = Mesh(shellSterileAdapter, c="red", alpha=0.4)
PSM2_SterileAdapter = Mesh(shellSterileAdapter, c="green", alpha=0.4)

FOV_depth = 0.2
FOV_width = 0.1
ECM_FOV = Cone(pos=(0,0,FOV_depth/2),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="cyan",alpha = 0.6)
ECM_FOV = ECM_FOV.pos(0,0,-FOV_depth/2)


mesh_algo_radius = 1
cloudECM = pv.PolyData(ECM_xyz)
#cloudECM.plot()
volumeECM = cloudECM.delaunay_3d(alpha = mesh_algo_radius)
shellECM = volumeECM.extract_geometry()
#shellECM.plot()

ECM_arm = Mesh(shellECM, c="blue", alpha=0.4)

write(ECM_FOV, CurrentDirectory + '/dVRK_Meshes/ECM_FOV.stl',binary=True)
write(Mesh(shellRWS), CurrentDirectory + '/dVRK_Meshes/PSM_EE_RWS.stl',binary=True)
write(Mesh(shellSterileAdapter), CurrentDirectory + '/dVRK_Meshes/PSM_SterileAdapter_Sweep.stl',binary=True)
write(Mesh(shellECM), CurrentDirectory + '/dVRK_Meshes/ECM_Arm_Sweep.stl',binary=True)

### Transformation of Volumes according to RCM locations ###
offset = 0.1 #for debugging RCM offset
rot_angle = 20 #for debugging RCM rotation deg

RCM1xyz = [offset,0,0] #[x,y,z] 
RCM2xyz = [-offset,0,0]
RCM_ECMxyz = [offset,0,0]

PSM1_RWS.pos(RCM1xyz).rotate(rot_angle,axis=(0,1,0),point=(RCM1xyz),rad=False)
PSM2_RWS.pos(RCM2xyz).rotate(-rot_angle,axis=(0,1,0),point=(RCM2xyz),rad=False)
PSM1_SterileAdapter.pos(RCM1xyz).rotate(rot_angle,axis=(0,1,0),point=(RCM1xyz),rad=False)
PSM2_SterileAdapter.pos(RCM2xyz).rotate(-rot_angle,axis=(0,1,0),point=(RCM2xyz),rad=False)
#ECM_FOV.pos(RCM_ECMxyz).rotate(20,axis=(0,1,0),point=(RCM_ECMxyz),rad=False)
ECM_arm.pos(RCM_ECMxyz).rotate(20,axis=(0,1,0),point=(RCM_ECMxyz),rad=False)

settings.use_depth_peeling = True
plt = Plotter(shape=(1,3), interactive=False, axes=3)
plt.at(0).show(PSM1_RWS, PSM1_SterileAdapter, PSM2_RWS, PSM2_SterileAdapter, ECM_FOV, ECM_arm, "environment", axes = True)
plt.interactive().close()

### intersect calculation ###
"""Opt1 = PSM1_RWS.boolean("intersect", ECM_FOV).c('magenta')
Opt1 = Opt1.boolean("intersect", PSM2_RWS).c('magenta')
print(Opt1.volume()*10**6 , "[cc]")
plt.at(1).show(Opt1, "intersect volume: %4.2f[cc]" % (Opt1.volume()*10**6) , resetcam=False)

ExtSweep1 = PSM1_SterileAdapter.boolean("intersect", ECM_arm).c('cyan')
ExtSweep2 = PSM2_SterileAdapter.boolean("intersect", ECM_arm).c('cyan')
Opt2 = ExtSweep1.boolean("intersect", PSM2_SterileAdapter).c('cyan')

if Opt2.volume()*10**6 > 0:
    intersect = Opt2.volume()*10**6
    print(intersect, "[cc]")
    plt.at(2).show(Opt2, "intersect volume: %4.2f[cc]" %intersect , resetcam=False)
else:
    intersect = (ExtSweep1.volume()+ExtSweep2.volume())*10**6 
    print(intersect, "[cc]")
    plt.at(2).show(ExtSweep1,ExtSweep2, "intersect volume: %4.2f[cc]" %intersect , resetcam=False)"""
#plt.interactive().close()