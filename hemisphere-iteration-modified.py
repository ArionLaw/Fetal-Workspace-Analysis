from math import sin, cos, pi
import math
import numpy as np
from vedo import*
import sys,time

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
plt = Plotter(shape=(1,1), interactive=False, axes=3)
port_locs = Mesh(zones).ps(10) # meshing and displaying points on hemisphere surface
# Generating hemisphere surface for boolean intersection
#uterus_boundary = Mesh(fit_sphere(zones))
Ext_Uterus_Simp = Sphere(pos=(0,0,0),r=R,alpha=0.4)  # generate full sphere
Ext_Uterus_Simp = Ext_Uterus_Simp.cut_with_plane(normal=(0,0,1)) # remove bottom of sphere
Ext_Uterus_Simp = Ext_Uterus_Simp.fill_holes(size=30) # fills bottom to make enclosed hemisphere
FOV_depth = 0.15
FOV_width = 0.05
ECM_FOV = Cone(pos=(0,0,0),r = FOV_width , height = FOV_depth, res=64, axis=(0,0,1),c="blue",alpha = 0.6)
show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1, azimuth = 30 , elevation = -30, roll = -30)   
"""ECM_FOV.pos(port[0,0] - [0,0,FOV_depth/2]).rotate(float(rot_ang[0,0]),axis=rot_axis[0,0], point=port[0,0],rad=True)
show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1)
plt.interactive().close()
ECM_FOV.rotate(-float(rot_ang[0,0]),axis=rot_axis[0,0], point=port[0,0],rad=True)

ECM_FOV.pos(port[0,3] - [0,0,FOV_depth/2]).rotate(float(rot_ang[0,3]),axis=rot_axis[0,3], point=port[0,3],rad=True)
show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1)
plt.interactive().close()
ECM_FOV.rotate(-float(rot_ang[0,3]),axis=rot_axis[0,3], point=port[0,3],rad=True)

ECM_FOV.pos(port[3,3] - [0,0,FOV_depth/2])
show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1)
plt.interactive().close()

ECM_FOV.pos(port[0,2] - [0,0,FOV_depth/2]).rotate(float(rot_ang[0,2]),axis=rot_axis[0,2], point=port[0,2],rad=True)
show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1)
plt.interactive().close()
ECM_FOV.rotate(-float(rot_ang[0,2]),axis=rot_axis[0,2], point=port[0,2],rad=True)"""
#print(port)
#print(rot_ang)
#print(rot_axis)
for a in range(steps):
    for b in range(steps):
        time.sleep(.1) #slowdown for visualization 
        if math.isnan(rot_axis[a,b,0]) == True or math.isnan(rot_axis[a,b,1]) == True or math.isnan(rot_axis[a,b,2]) == True:
            ECM_FOV.pos(port[a,b] - [0,0,FOV_depth/2])
            show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1)            
        else:
            ECM_FOV.pos(port[a,b] - [0,0,FOV_depth/2]).rotate(float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)
            show(port_locs, Ext_Uterus_Simp, ECM_FOV, __doc__, axes=1)
            ECM_FOV.rotate(-float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)
ECM_FOV.rotate(float(rot_ang[a,b]),axis=rot_axis[a,b], point=port[a,b],rad=True)            
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