from math import sin, cos, pi
import math
import numpy as np
from vedo import Plotter, Mesh, Sphere, Cone

def MeshTransformtoPort(meshInternal,meshExternal,a,b):
    if math.isnan(rot_axis[a,b,0]) == True: #all values NaN if no rotation
        meshInternal.pos(port[a,b] - [0,0,FOV_depth/2])
        meshExternal.pos(port[a,b] + [0,0,FOV_depth/2])    
    else:
        meshInternal.pos(port[a,b] - [0,0,FOV_depth/2]).rotate(float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)
        meshExternal.pos(port[a,b] + [0,0,FOV_depth/2]).rotate(float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)

def MeshReset(meshInternal,meshExternal,a,b):
    if math.isnan(rot_axis[a,b,0]) == True: #all values NaN if no rotation
        meshInternal.pos(port[a,b] - [0,0,FOV_depth/2])
        meshExternal.pos(port[a,b] + [0,0,FOV_depth/2])    
    else:
        meshInternal.rotate(-float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)
        meshExternal.rotate(-float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)

def getInternalIntersect(ECM_FOV, PSM1_RWS , PSM2_RWS):
    ECM_PSM1 = ECM_FOV.boolean("intersect", PSM1_RWS).c('magenta')
    Opt1 = ECM_PSM1.boolean("intersect", PSM2_RWS).c('magenta')
    intersect = Opt1.volume()*10**6
    print("internal intersect volume: %4.2f[cc]" %intersect)
    plt.at(1).show(Opt1, "internal intersect volume: %4.2f[cc]" %intersect , resetcam=False)
    return intersect

def getExternalIntersect(ECM_arm, PSM1_SterileAdapter , PSM2_SterileAdapter):
    ExtSweep1 = ECM_arm.boolean("intersect", PSM1_SterileAdapter).c('cyan')
    ExtSweep2 = ECM_arm.boolean("intersect", PSM2_SterileAdapter).c('cyan')
    Opt2 = ExtSweep1.boolean("intersect", PSM2_SterileAdapter).c('cyan')

    if Opt2.volume()*10**6 > 0:
        intersect = Opt2.volume()*10**6
        print("external intersect volume: %4.2f[cc]" %intersect)
        plt.at(2).show(Opt2, "external intersect volume: %4.2f[cc]" %intersect , resetcam=False)
        return intersect
    else:
        intersect = (ExtSweep1.volume()+ExtSweep2.volume())*10**6 
        print("external intersect volume: %4.2f[cc]" %intersect)
        plt.at(2).show(ExtSweep1,ExtSweep2, "external intersect volume: %4.2f[cc]" %intersect, resetcam=False)
        return intersect

R = 0.15 # radius (m)
i = 0
min = -30*pi/180    # -90 deg to rad
max = 30*pi/180     # 90 deg to rad
step = 10*pi/180     # 5 deg step size to rad
steps = int((max-min)/step)+1

alpha = np.linspace(min,max,steps)
beta = np.linspace(min,max,steps)
zones = np.ndarray(shape=(steps**2,3)) # array to store surface points
port = np.ndarray(shape=(steps,steps,3)) # matrix to store surface points
rot_ang = np.ndarray(shape=(steps,steps,1)) # matrix to store normal angle
rot_axis = np.ndarray(shape=(steps,steps,3)) # matrix to store normal axis of rotation

# generating hemisphere using spherical coords
for a in range(steps):
    for b in range(steps):
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

plt = Plotter(shape=(1,3), interactive=False, axes=1)

port_locs = Mesh(zones).ps(10) # meshing and displaying points on hemisphere surface
Ext_Uterus_Simp = Sphere(pos=(0,0,0),r=R,c="red",alpha=0.4)  # generate full sphere
Ext_Uterus_Simp = Ext_Uterus_Simp.cut_with_plane(normal=(0,0,1)) # remove bottom of sphere
Ext_Uterus_Simp = Ext_Uterus_Simp.fill_holes(size=30) # fills bottom to make enclosed hemisphere

FOV_depth = 0.15
FOV_width = 0.08
ECM_FOV = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="cyan",alpha = 0.6)
PSM1_EE = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="magenta",alpha = 0.6)
PSM2_EE = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="green",alpha = 0.6)
ECM_sweep = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="cyan",alpha = 0.6)
PSM1_sweep = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="magenta",alpha = 0.6)
PSM2_sweep = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="green",alpha = 0.6)
ECM_sweep.rotate(180,axis=[0,1,0], point=[0,0,0],rad=False)
PSM1_sweep.rotate(180,axis=[0,1,0], point=[0,0,0],rad=False)
PSM2_sweep.rotate(180,axis=[0,1,0], point=[0,0,0],rad=False)



partitionLS_PSM1 = range(steps//2)
partitionRS_PSM2 = range(steps//2+1 , steps)
partitionECM = [steps//2]

sampleSize = steps*len(partitionLS_PSM1)*steps*len(partitionRS_PSM2)*steps
Optimal_Intersect_Data = [[]] #store a1,b1,a2,b2,a3,b3, internal intersect , external intersect

"""a1 = 1
b1 = 1
a2 = 2
b2 = 1
a3 = 3
b3 = 1
MeshTransformtoPort(ECM_FOV,ECM_sweep,a1,b1)
MeshTransformtoPort(PSM1_EE,PSM1_sweep,a2,b2)
MeshTransformtoPort(PSM2_EE,PSM2_sweep,a3,b3)
plt.at(0).show(port_locs, Ext_Uterus_Simp, ECM_FOV, PSM1_EE , PSM2_EE , ECM_sweep, PSM1_sweep , PSM2_sweep, __doc__, axes=1, camera = {'pos':(0.5,-0.5,1), 'focal_point':(0,0,0), 'viewup':(0,0,1)})

print("estimated run time: %5i / %5i" %(i , sampleSize))
internal = getInternalIntersect(ECM_FOV,PSM1_EE,PSM2_EE)
external = getExternalIntersect(ECM_sweep,PSM1_sweep,PSM2_sweep)
Optimal_Intersect_Data.append([a,b,a2,b2,a3,b3,internal,external])

MeshReset(PSM2_EE,PSM2_sweep,a3,b3)
MeshReset(PSM1_EE,PSM1_sweep,a2,b2)
MeshReset(ECM_FOV,ECM_sweep,a1,b1)"""

"""a1 = 1
b1 = 0
a2 = 2
b2 = 1
a3 = 3
b3 = 2
MeshTransformtoPort(ECM_FOV,ECM_sweep,a1,b1)
MeshTransformtoPort(PSM1_EE,PSM1_sweep,a2,b2)
MeshTransformtoPort(PSM2_EE,PSM2_sweep,a3,b3)
plt.at(0).show(port_locs, Ext_Uterus_Simp, ECM_FOV, PSM1_EE , PSM2_EE , ECM_sweep, PSM1_sweep , PSM2_sweep, __doc__, axes=1, camera = {'pos':(0.5,-0.5,1), 'focal_point':(0,0,0), 'viewup':(0,0,1)})

print("estimated run time: %5i / %5i" %(i , sampleSize))
internal = getInternalIntersect(ECM_FOV,PSM1_EE,PSM2_EE)
external = getExternalIntersect(ECM_sweep,PSM1_sweep,PSM2_sweep)
Optimal_Intersect_Data.append([a,b,a2,b2,a3,b3,internal,external])

plt.interactive().close()
MeshReset(PSM2_EE,PSM2_sweep,a3,b3)
MeshReset(PSM1_EE,PSM1_sweep,a2,b2)
MeshReset(ECM_FOV,ECM_sweep,a1,b1)"""

i=0
for a1 in partitionECM:
    for b1 in range(steps):
        MeshTransformtoPort(ECM_FOV,ECM_sweep,a1,b1)       
        
        for a2 in partitionLS_PSM1:
            for b2 in range(steps):
                MeshTransformtoPort(PSM1_EE,PSM1_sweep,a2,b2)
                
                for a3 in partitionRS_PSM2:
                    for b3 in range(steps):
                        MeshTransformtoPort(PSM2_EE,PSM2_sweep,a3,b3)
                        plt.at(0).show(port_locs, Ext_Uterus_Simp, ECM_FOV, PSM1_EE , PSM2_EE , ECM_sweep, PSM1_sweep , PSM2_sweep, __doc__, axes=1, camera = {'pos':(0.5,-0.5,1), 'focal_point':(0,0,0), 'viewup':(0,0,1)})
                        #getInternalIntersect(ECM_FOV,PSM1_EE,PSM2_EE)
                        #getExternalIntersect(ECM_sweep,PSM1_sweep,PSM2_sweep)
                        """print("estimated run time: %5i / %5i" %(i , sampleSize))
                        internal = getInternalIntersect(ECM_FOV,PSM1_EE,PSM2_EE)
                        external = getExternalIntersect(ECM_sweep,PSM1_sweep,PSM2_sweep)
                        print("internal intersect: %4.2f [cc]" %internal)
                        print("external intersect: %4.2f [cc]" %external)
                        Optimal_Intersect_Data.append([a,b,a2,b2,a3,b3,internal,external])"""
                        #plt.interactive().close()
                        i+=1
                        
                        MeshReset(PSM2_EE,PSM2_sweep,a3,b3)
                MeshReset(PSM1_EE,PSM1_sweep,a2,b2)
        MeshReset(ECM_FOV,ECM_sweep,a1,b1)
                                    
ECM_FOV.rotate(float(rot_ang[a1,b1]),axis=rot_axis[a1,b1], point=port[a1,b1],rad=True)
PSM1_EE.rotate(float(rot_ang[a2,b2]),axis=rot_axis[a2,b2], point=port[a2,b2],rad=True)
PSM2_EE.rotate(float(rot_ang[a3,b3]),axis=rot_axis[a3,b3], point=port[a3,b3],rad=True)
ECM_sweep.rotate(float(rot_ang[a1,b1]),axis=rot_axis[a1,b1], point=port[a1,b1],rad=True)
PSM1_sweep.rotate(float(rot_ang[a2,b2]),axis=rot_axis[a2,b2], point=port[a2,b2],rad=True)
PSM2_sweep.rotate(float(rot_ang[a3,b3]),axis=rot_axis[a3,b3], point=port[a3,b3],rad=True)
plt.interactive().close()


# Everything below this point is pseudocode - comment out if you want to see the above code run
"""# Import 3 meshes from the 3 arms
# psm1_mesh = //mesh object #1: PSM1 ROM mesh
# psm2_mesh = //mesh object #2: PSM2 ROM mesh
# ecm_mesh = ///mesh object #3: ECM ROM mesh

# Iterating through surface points x3
    # may be more efficient to drop the prior sphere construction and include that within these nested loops
    # as there is no reason to construct if we repeat that iteration anyways...
for a1 in np.arange(min,max,step):
    for b1 in np.arange(min,max,step):
        # align PSM1 mesh RCM with this point:
        # RCM1 at (a1,b1)
        # PSM1 normal aligned with normal1 = (sin(a1)*cos(b1),sin(a1)*sin(b1),cos(b1)) (normal to sphere at surface point)

        for a2 in np.arange(min, max, step):
            for b2 in np.arange(min, max, step):
                # align PSM2 mesh RCM with (a2,b2), normal2
                # if (a2,b2) == (a1, b1), break (move on to next b2 so both PSMs are not at same point)

                for a3 in np.arange(min, max, step):
                    for b3 in np.arange(min, max, step):
                        # align ECM mesh RCM with (a3,b3), normal3
                        # if (a3,b3) == (a1, b1) || (a3,b3) == (a2,b2), break (skip to next pt)

                        # INNER WORKINGS:

                        #slice PSM1_mesh by uterus_boundary, into:
                            #PSM1_outer
                            #PSM1_inner
                        #slice PSM2_mesh by uterus_boundary, into:
                            # PSM2_outer
                            # PSM2_inner
                        #slice ECM_mesh by uterus_boundary, into:
                            # ECM_outer
                            # ECM_inner

                        #workspace = PSM1_inner.boolean("intersect", PSM2_inner)
                        #visible_workspace = workspace.boolean("intersect", ECM_inner)
                        #c1 = PSM1_outer.boolean("intersect", PSM2_outer)
                        #c2 = PSM1_outer.boolean("intersect", ECM_outer)
                        #c3 = PSM2_outer.boolean("intersect", ECM_outer)
                        # collisions = [c1,c2,c3]

                        #store all these values:
                        # [a1,b1,a2,b2,a3,b3,visible_workspace,collisions[]]
                        # New script to assign optimization score to each array element and sort to find 'best' port placements
                            # Compare visible_workspace vol : total collisions vol
                            # Calculate port distances between each other, from center uterus point (as informed by Tim)
# END of 6x nested loops"""