from math import sin, cos, pi
import numpy as np
from vedo import Mesh, show, Volume, Sphere, fit_sphere

length = 36*36  # no. of steps ^2
zones = np.ndarray(shape=(length,3)) # array to store surface points

R = 15 # radius (cm)
i = 0
min = -90*pi/180    # -90 deg to rad
max = 90*pi/180     # 90 deg to rad
step = 5*pi/180     # 5 deg step size to rad

# generating hemisphere using spherical coords
for a in np.arange(min,max,step):
    for b in np.arange(min,max,step):
        x = R*sin(a)*cos(b)
        y = R*sin(a)*sin(b)
        z = R*cos(a)

        zones[i] = (x,y,z)
        i += 1

insertions = Mesh(zones).ps(10) # meshing and displaying points on hemisphere surface

# Generating hemisphere surface for boolean intersection
uterus_boundary = Mesh(fit_sphere(zones))
sph = Sphere(pos=(0,0,0),r=15)  # generate full sphere
cut_sphere = sph.cut_with_plane(normal=(0,0,1)) # remove bottom of sphere
internal = cut_sphere.fill_holes(size=30) # fills bottom to make enclosed hemisphere
# there is probably a cleaner way to fill the bottom
# need some method of filling the sphere for the purpose of future bool intersections

show(insertions, internal, __doc__, viewup='z', axes=1) # temp: visualizing constructed uterus

# Everything below this point is pseudocode - comment out if you want to see the above code run

# Import 3 meshes from the 3 arms
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
# END of 6x nested loops