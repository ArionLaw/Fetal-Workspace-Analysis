from math import sin, cos, pi
import math
import numpy as np
from vedo import Plotter, Mesh, Sphere, Ellipsoid, Line
import time
import os
import sys

CurrentDirectory = sys.path[0]

def Magnitude(vector):
    return math.sqrt(sum(pow(element, 2) for element in vector))

def getPortToLesionData(a1,b1,a2,b2,a3,b3):
    #surface normals at port locations and fetus lesion 
    LesionNorm = LesionNormalVector.pos() - Lesion.pos()
    ECMPortNorm = port_norm[a1,b1]
    PSM1PortNorm = port_norm[a2,b2]
    PSM2PortNorm = port_norm[a3,b3]
    #vectors from ports to lesion
    ECMtoLesion = port[a1,b1] - Lesion.pos()
    PSM1toLesion = port[a2,b2] - Lesion.pos()
    PSM2toLesion = port[a3,b3] - Lesion.pos()
    #distance from ports to lesion
    ECMtoLesionDistance = Magnitude(ECMtoLesion)
    PSM1toLesionDistance = Magnitude(PSM1toLesion)
    PSM2toLesionDistance = Magnitude(PSM2toLesion)
    #difference in angle from vectors to lesion normal
    ECMtoLesionNormAngle = np.arccos(np.dot(ECMtoLesion,LesionNorm)/ECMtoLesionDistance/Magnitude(LesionNorm))
    PSM1toLesionNormAngle = np.arccos(np.dot(PSM1toLesion,LesionNorm)/PSM1toLesionDistance/Magnitude(LesionNorm))
    PSM2toLesionNormAngle = np.arccos(np.dot(PSM2toLesion,LesionNorm)/PSM2toLesionDistance/Magnitude(LesionNorm))
    #difference in angle from vectors to port normal
    ECMPortTwistAngle = np.arccos(np.dot(ECMtoLesion,ECMPortNorm)/ECMtoLesionDistance/Magnitude(ECMPortNorm))
    PSM1PortTwistAngle = np.arccos(np.dot(PSM1toLesion,PSM1PortNorm)/PSM1toLesionDistance/Magnitude(PSM1PortNorm))
    PSM2PortTwistAngle = np.arccos(np.dot(PSM2toLesion,PSM2PortNorm)/PSM2toLesionDistance/Magnitude(PSM2PortNorm))

    return [ECMtoLesionDistance, PSM1toLesionDistance, PSM2toLesionDistance, ECMtoLesionNormAngle, PSM1toLesionNormAngle, PSM2toLesionNormAngle, ECMPortTwistAngle, PSM1PortTwistAngle, PSM2PortTwistAngle]

def MeshTransformtoPort(mesh,a,b):
    if math.isnan(rot_axis[a,b,0]) == True: #all values NaN if no rotation
        mesh.pos(port[a,b]) 
    else:
        mesh.pos(port[a,b]).rotate(float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)
    
    PortNorm = port_norm[a,b]
    PortToLesion = port[a,b] - Lesion.pos()
    PortTwistAngle = np.arccos(np.dot(PortToLesion,PortNorm)/Magnitude(PortToLesion)/Magnitude(PortNorm))
    PortRotAxis = np.cross(PortToLesion,PortNorm)/Magnitude(PortToLesion)/Magnitude(PortNorm)
    
    mesh.pos(port[a,b]).rotate(float(PortTwistAngle),axis=-PortRotAxis, point=port[a,b],rad=True)

def MeshReset(mesh,a,b):
    PortNorm = port_norm[a,b]
    PortToLesion = port[a,b] - Lesion.pos()
    PortTwistAngle = np.arccos(np.dot(PortToLesion,PortNorm)/Magnitude(PortToLesion)/Magnitude(PortNorm))
    PortRotAxis = np.cross(PortToLesion,PortNorm)/Magnitude(PortToLesion)/Magnitude(PortNorm)
    
    mesh.pos(port[a,b]).rotate(-float(PortTwistAngle),axis=-PortRotAxis, point=port[a,b],rad=True)
    
    if math.isnan(rot_axis[a,b,0]) == True: #all values NaN if no rotation
        mesh.pos(port[a,b]) 
    else:
        mesh.pos(port[a,b]).rotate(-float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)    

def getInternalIntersect(ECM_FOV, PSM1_RWS , PSM2_RWS):
    #plt.at(1).clear(deep=True)
    ECM_U_PSM1 = ECM_FOV.boolean("intersect", PSM1_RWS).c('magenta')
    ECM_U_PSM1_U_PSM2 = ECM_U_PSM1.boolean("intersect", PSM2_RWS).c('magenta')
    intersect = ECM_U_PSM1_U_PSM2.volume()*10**6
    print("internal intersect volume: %4.2f[cc]" %intersect)
    #plt.at(1).show(ECM_U_PSM1_U_PSM2, "internal intersect volume: %4.2f[cc]" %intersect , interactive=False, axes=1, resetcam=False)
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

### optimization domain ###
xmin = -15*pi/180
xmax = 15*pi/180
ymin = -30*pi/180
ymax = 20*pi/180
fmin = 0
fmax = 30*pi/180
step = 5*pi/180     # 5 deg step size to rad
xsteps = int((xmax-xmin)/step)+1
ysteps = int((ymax-ymin)/step)+1
fsteps = int((fmax-fmin)/step)+1
alpha = np.linspace(xmin,xmax,xsteps)
beta = np.linspace(ymin,ymax,ysteps)
theta = np.linspace(fmin,fmax,fsteps)
zones = np.ndarray(shape=(xsteps*ysteps,3)) # array to store surface points
port = np.ndarray(shape=(xsteps,ysteps,3)) # matrix to store surface points
port_norm = np.ndarray(shape=(xsteps,ysteps,3)) # matrix to store surface points
rot_ang = np.ndarray(shape=(xsteps,ysteps,1)) # matrix to store normal angle
rot_axis = np.ndarray(shape=(xsteps,ysteps,3)) # matrix to store normal axis of rotation

### generating hemisphere using spherical coords ###
R = 0.13 # uterus radius (m)
i = 0
for a in range(xsteps):
    for b in range(ysteps):
        x = R*sin(alpha[a])#*cos(beta[b])
        y = R*cos(alpha[a])*sin(beta[b])
        z = R*cos(alpha[a])*cos(beta[b])
        zones[i] = [x,y,z]
        
        v1_u = [x,y,z]/np.linalg.norm([x,y,z])
        v2_u = [x,y,0]/np.linalg.norm([x,y,0])
        v3_u = [0,0,1]

        port[a,b] = zones[i]
        port_norm[a,b] = v1_u
        rot_ang[a,b] = np.arccos(np.dot(v1_u, v3_u)/Magnitude(v1_u)/Magnitude(v3_u))
        rot_axis[a,b] = np.cross(v1_u,v2_u)/Magnitude(v1_u)/Magnitude(v2_u)
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
LesionNormalVector = Line([0,0.045,0.046],[0,0.045,0.086],c="red",lw=2).pattern('- -',repeats = 10)

ECM_FOV = Mesh(CurrentDirectory + '/dVRK_Meshes/ECM_FOV.stl',c = "cyan", alpha=0.6)
PSM1_EE = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_EE_RWS.stl',c = "magneta", alpha=0.6)
PSM2_EE = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_EE_RWS.stl',c = "green", alpha=0.6)
#ECM_sweep = Mesh(CurrentDirectory + '/dVRK_Meshes/ECM_Arm_Sweep.stl',c = "cyan", alpha=0.2)
#PSM1_sweep = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_SterileAdapter_Sweep.stl',c = "magneta", alpha=0.2)
#PSM2_sweep = Mesh(CurrentDirectory + '/dVRK_Meshes/PSM_SterileAdapter_Sweep.stl',c = "green", alpha=0.2)

partitionLS_PSM1 = range(0,xsteps//2)
partitionRS_PSM2 = range(xsteps,xsteps//2+1)
partitionECM = [xsteps//2]
TotalSamples = fsteps*ysteps*len(partitionLS_PSM1)
OptimizationData = [[]] 

#"""
offset = 3
spinalLocation = 0
FetalPivotAngle = 30
Fetus.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=False)
Lesion.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=False)
LesionNormalVector.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=False)
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

internal = getInternalIntersect(ECM_FOV = ECM_FOV , PSM1_RWS = PSM1_EE , PSM2_RWS=PSM2_EE)
OptimizationData.append([a1,b1,a2,b2,a3,b3,internal] + getPortToLesionData(a1,b1,a2,b2,a3,b3))
plt.at(0).show(port_locs, Ext_Uterus_Simp, Uterus, Fetus, Lesion, LesionNormalVector, ECM_FOV, PSM1_EE , PSM2_EE, __doc__, axes=1, camera = {'pos':(0.3,-0.6,0.6), 'focal_point':(0,0,0), 'viewup':(0,0,1)})

print(OptimizationData)
plt.interactive().close()

MeshReset(ECM_FOV,a1,b1)
MeshReset(PSM1_EE,a2,b2)
MeshReset(PSM2_EE,a3,b3)
#MeshReset(ECM_sweep,a1,b1)
#MeshReset(PSM1_sweep,a2,b2)
#MeshReset(PSM2_sweep,a3,b3)
#"""

i=0
for t in range(fsteps):
    FetalPivotAngle = theta[t]
    Fetus.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
    Lesion.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
    LesionNormalVector.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
    for a1 in partitionECM:
        for b1 in range(ysteps):
            MeshTransformtoPort(ECM_FOV,a1,b1)

            for a2 in partitionLS_PSM1:
                b2 = b1
                a3 = partitionECM[0] + (partitionECM[0]-a2)
                b3 = b1
                MeshTransformtoPort(PSM1_EE,a2,b2)  
                MeshTransformtoPort(PSM2_EE,a3,b3)

                #time.sleep(0.2)
                plt.at(0).show(port_locs, Ext_Uterus_Simp, Uterus, Fetus, Lesion, LesionNormalVector, ECM_FOV, PSM1_EE , PSM2_EE, __doc__, axes=1, camera = {'pos':(0.3,-0.6,0.6), 'focal_point':(0,0,0), 'viewup':(0,0,1)})
                print("estimated run time: %5i / %5i" %(i , TotalSamples))
                #internal = getInternalIntersect(ECM_FOV = ECM_FOV , PSM1_RWS = PSM1_EE , PSM2_RWS=PSM2_EE)
                #OptimizationData.append([t,a1,b1,a2,b2,a3,b3,internal] + getPortToLesionData(a1,b1,a2,b2,a3,b3))
                OptimizationData.append([t,a1,b1,a2,b2,a3,b3] + getPortToLesionData(a1,b1,a2,b2,a3,b3))
                #plt.interactive().close()
                i+=1

                MeshReset(PSM2_EE,a3,b3)
                MeshReset(PSM1_EE,a2,b2)
            MeshReset(ECM_FOV,a1,b1)
    Fetus.rotate(-FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
    Lesion.rotate(-FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
    LesionNormalVector.rotate(-FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
 
Fetus.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
Lesion.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
LesionNormalVector.rotate(FetalPivotAngle,axis=(1,0,0),point=FetalHeadPivot,rad=True)
MeshTransformtoPort(ECM_FOV,a1,b1)
MeshTransformtoPort(PSM1_EE,a2,b2)  
MeshTransformtoPort(PSM2_EE,a3,b3)

OptimizationData.pop(0)
print(OptimizationData)
plt.interactive().close()