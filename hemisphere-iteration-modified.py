from math import sin, cos, pi
import math
import numpy as np
from vedo import Plotter, Mesh, Sphere, Ellipsoid
import time
import os
import sys

CurrentDirectory = sys.path[0]

def MeshTransformtoPort(mesh,a,b):
    if math.isnan(rot_axis[a,b,0]) == True: #all values NaN if no rotation
        mesh.pos(port[a,b]) 
    else:
        mesh.pos(port[a,b]).rotate(float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)

def MeshReset(mesh,a,b):
    if math.isnan(rot_axis[a,b,0]) == True: #all values NaN if no rotation
        mesh.pos(port[a,b]) 
    else:
        mesh.pos(port[a,b]).rotate(-float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)

def getInternalIntersect(ECM_FOV, PSM1_RWS , PSM2_RWS):
    ECM_PSM1 = ECM_FOV.boolean("intersect", PSM1_RWS).c('magenta')
    Opt1 = ECM_PSM1.boolean("intersect", PSM2_RWS).c('magenta')
    intersect = Opt1.volume()*10**6
    print("internal intersect volume: %4.2f[cc]" %intersect)
    #plt.at(1).show(Opt1, "internal intersect volume: %4.2f[cc]" %intersect , resetcam=False)
    return intersect

"""
def getExternalIntersect(ECM_arm, PSM1_SterileAdapter , PSM2_SterileAdapter):
    ExtSweep1 = ECM_arm.boolean("intersect", PSM1_SterileAdapter).c('cyan')
    ExtSweep2 = ECM_arm.boolean("intersect", PSM2_SterileAdapter).c('cyan')
    Opt2 = ExtSweep1.boolean("intersect", PSM2_SterileAdapter).c('cyan')

    #fix collision defn 

    if Opt2.volume()*10**6 > 0:
        intersect = Opt2.volume()*10**6
        print("external intersect volume: %4.2f[cc]" %intersect)
        plt.at(2).show(Opt2, "external intersect volume: %4.2f[cc]" %intersect , resetcam=False)
    else:
        intersect = (ExtSweep1.volume()+ExtSweep2.volume())*10**6 
        print("external intersect volume: %4.2f[cc]" %intersect)
        plt.at(2).show(ExtSweep1,ExtSweep2, "external intersect volume: %4.2f[cc]" %intersect, resetcam=False)
    
    return intersect
"""
xmin = -15*pi/180
xmax = 15*pi/180
ymin = -45*pi/180
ymax = 45*pi/180
step = 5*pi/180     # 5 deg step size to rad
xsteps = int((xmax-xmin)/step)+1
ysteps = int((ymax-ymin)/step)+1
alpha = np.linspace(xmin,xmax,xsteps)
beta = np.linspace(ymin,ymax,ysteps)
zones = np.ndarray(shape=(xsteps*ysteps,3)) # array to store surface points
port = np.ndarray(shape=(xsteps,ysteps,3)) # matrix to store surface points
rot_ang = np.ndarray(shape=(xsteps,ysteps,1)) # matrix to store normal angle
rot_axis = np.ndarray(shape=(xsteps,ysteps,3)) # matrix to store normal axis of rotation

# generating hemisphere using spherical coords
R = 0.13 # uterus radius (m)
i = 0
for a in range(xsteps):
    for b in range(ysteps):
        x = R*sin(alpha[a])#*cos(beta[b])
        y = R*cos(alpha[a])*sin(beta[b])
        z = R*cos(alpha[a])*cos(beta[b])
        zones[i] = [x,y,z]
        port[a,b] = zones[i]
        v1_u = [x,y,z]/np.linalg.norm([x,y,z])
        v2_u = [x,y,0]/np.linalg.norm([x,y,0])
        v3_u = [0,0,1]
        rot_ang[a,b] = np.arccos(np.clip(np.dot(v1_u, v3_u), -1.0, 1.0))
        rot_axis[a,b] = np.cross(v1_u,v2_u)/np.linalg.norm(np.cross(v1_u,v2_u))
        i += 1

plt = Plotter(shape=(1,1), interactive=False, axes=1)

port_locs = Mesh(zones).ps(10) # meshing and displaying points on hemisphere surface
Ext_Uterus_Simp = Sphere(pos=(0,0,0),r=R,c="red",alpha=0.4)  # generate full sphere
Ext_Uterus_Simp = Ext_Uterus_Simp.cut_with_plane(normal=(0,0,1)) # remove bottom of sphere
Ext_Uterus_Simp = Ext_Uterus_Simp.fill_holes(size=30) # fills bottom to make enclosed hemisphere

Uterus = Mesh(CurrentDirectory + '/Fetal_Meshes/UterusPhantom.stl',c = "pink", alpha=0.4)
Uterus.rotate(-90,axis=(1,0,0),point=(0,0,0),rad=False).scale(s=0.0013,reset=True).shift(dx=0,dy=0.048,dz=0.048)
Fetus = Mesh(CurrentDirectory + '/Fetal_Meshes/FetalPhantom.stl',c = "yellow", alpha=0.6)
Fetus.rotate(-90,axis=(1,0,0),point=(0,0,0),rad=False).rotate(180,axis=(0,0,1),point=(0,0,0),rad=False).scale(s=0.001,reset=True).shift(dx=0,dy=0.01,dz=0.06)
FetalHeadPivot = [0,-0.05,0]
Lesion = Ellipsoid(pos=(0,0.045,0.046),axis1=(0.018,0,0),axis2=(0,0.028,0),axis3=(0,0,0.006),res=24,c="red",alpha=0.8)

ECM_FOV = Mesh(CurrentDirectory + '/dVRK_Meshes/ECM_FOV.stl',c = "cyan", alpha=0.6)
PSM1_EE = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_EE_RWS.stl',c = "magneta", alpha=0.6)
PSM2_EE = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_EE_RWS.stl',c = "green", alpha=0.6)
#ECM_sweep = Mesh(CurrentDirectory + '/dVRK_Meshes/ECM_Arm_Sweep.stl',c = "cyan", alpha=0.2)
#PSM1_sweep = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_SterileAdapter_Sweep.stl',c = "magneta", alpha=0.2)
#PSM2_sweep = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_SterileAdapter_Sweep.stl',c = "green", alpha=0.2)

partitionLS_PSM1 = range(0,xsteps//2)
partitionRS_PSM2 = range(xsteps,xsteps//2+1)
partitionECM = [xsteps//2]

TotalSamples = ysteps*len(partitionLS_PSM1)
Optimal_Intersect_Data = [[]] 

offset = 2
spinalLocation = 8
FetalPivotAngle = 20
Fetus.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=False)
Lesion.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=False)
a1 = partitionECM[0]
b1 = spinalLocation
a2 = partitionECM[0]-offset
b2 = spinalLocation
a3 = partitionECM[0]+offset
b3 = spinalLocation
MeshTransformtoPort(ECM_FOV,a1,b1)
MeshTransformtoPort(PSM1_EE,a2,b2)
MeshTransformtoPort(PSM2_EE,a3,b3)
#MeshTransformtoPort(ECM_sweep,a1,b1)
#MeshTransformtoPort(PSM1_sweep,a2,b2)
#MeshTransformtoPort(PSM2_sweep,a3,b3)

internal = getInternalIntersect(ECM_FOV,PSM1_EE,PSM2_EE)
#external = getExternalIntersect(ECM_sweep,PSM1_sweep,PSM2_sweep)
Optimal_Intersect_Data.append([a,b,a2,b2,a3,b3,internal])
plt.at(0).show(port_locs, Ext_Uterus_Simp, Uterus, Fetus, Lesion, ECM_FOV, PSM1_EE , PSM2_EE, __doc__, axes=1, camera = {'pos':(0.3,-0.6,0.6), 'focal_point':(0,0,0), 'viewup':(0,0,1)})
plt.interactive().close()

MeshReset(ECM_FOV,a1,b1)
MeshReset(PSM1_EE,a2,b2)
MeshReset(PSM2_EE,a3,b3)
#MeshReset(ECM_sweep,a1,b1)
#MeshReset(PSM1_sweep,a2,b2)
#MeshReset(PSM2_sweep,a3,b3)

i=0
for a1 in partitionECM:
    for b1 in range(ysteps):
        MeshTransformtoPort(ECM_FOV,a1,b1)

        for a2 in partitionLS_PSM1:
            b2 = b1
            a3 = partitionECM[0] + (partitionECM[0]-a2)
            b3 = b1
            MeshTransformtoPort(PSM1_EE,a2,b2)  
            MeshTransformtoPort(PSM2_EE,a3,b3)

            #time.sleep(0.1)
            plt.at(0).show(port_locs, Ext_Uterus_Simp, Uterus, Fetus, ECM_FOV, PSM1_EE , PSM2_EE, __doc__, axes=1, camera = {'pos':(0.3,-0.6,0.6), 'focal_point':(0,0,0), 'viewup':(0,0,1)})
            print("estimated run time: %5i / %5i" %(i , TotalSamples))
            internal = getInternalIntersect(ECM_FOV,PSM1_EE,PSM2_EE)
            Optimal_Intersect_Data.append([a,b,a2,b2,a3,b3,internal])
            #plt.interactive().close()
            i+=1

            MeshReset(PSM2_EE,a3,b3)
            MeshReset(PSM1_EE,a2,b2)
        MeshReset(ECM_FOV,a1,b1)
                    
ECM_FOV.rotate(float(rot_ang[a1,b1]),axis=rot_axis[a1,b1], point=port[a1,b1],rad=True)
PSM1_EE.rotate(float(rot_ang[a2,b2]),axis=rot_axis[a2,b2], point=port[a2,b2],rad=True)
PSM2_EE.rotate(float(rot_ang[a3,b3]),axis=rot_axis[a3,b3], point=port[a3,b3],rad=True)
plt.interactive().close()