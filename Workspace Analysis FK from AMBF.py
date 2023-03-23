from vedo import*
import numpy as np
from psmFK import*
import pyvista as pv


PI = np.pi

### For PSM ###
"""
q1 Outer Yaw
q2 Outer Pitch
q3 Insertion
q4 Wrist Roll
q5 Wrist Pitch
q6 Wrist Yaw
"""
# joint max limits: Outer Yaw, Outer Pitch, Insertion, Wrist Roll, Wrist Pitch, Wrist Yaw 
q_max = [90*PI/180 , 54*PI/180 , 0.240 , 180*PI/180 , 180*PI/180 , 180*PI/180 ] # in radians and meters
q_min = [-90*PI/180 , -54*PI/180 , 0.000 , -180*PI/180 , -180*PI/180 , -180*PI/180 ] # in radians and meters

# joint soft limits: Outer Yaw, Outer Pitch, Insertion, Wrist Roll, Wrist Pitch, Wrist Yaw 
qlim_u = [40*PI/180 , 30*PI/180 , 0.2 , 120*PI/180 , 90*PI/180 , 90*PI/180 ] # in radians and meters
qlim_l = [-40*PI/180 , -30*PI/180 , 0.02 , -120*PI/180 , -90*PI/180 , -90*PI/180 ] # in radians and meters

revolute1_divisions = 5
revolute2_divisions = 5
prismatic_divisions = round((qlim_u[2] - qlim_l[2])/0.005)+1
print(prismatic_divisions)

q1 = np.linspace(qlim_l[0],qlim_u[0],revolute1_divisions)
q2 = np.linspace(qlim_l[1],qlim_u[1],revolute2_divisions)
q3 = np.linspace(qlim_l[2],qlim_u[2],prismatic_divisions)

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
            #print("Outer Yaw [rads]: " , q1[h])
            #print("Outer Pitch [rads]: " , q2[k])
            #print("Insertion [m]: " , q3[j] )
            #print("End Effector Position[m]: ")
            #print( PSM_EE_pos , '\n')
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
    if PSM_EE_xyz[i][0] - PSM_EE_xyz[i-1][0] > 0.005:
        revolute1_divisions = revolute1_divisions + 2
        q1 = np.linspace(qlim_l[0],qlim_u[0],revolute1_divisions)
        revolute2_divisions = revolute2_divisions + 2
        q2 = np.linspace(qlim_l[1],qlim_u[1],revolute2_divisions)      

PSM_EE_xyz.pop(0)

mesh_algo_radius = 1

cloud = pv.PolyData(PSM_EE_xyz)
#cloud.plot()
volume = cloud.delaunay_3d(alpha = mesh_algo_radius)
shell = volume.extract_geometry()
#shell.plot()
PSM1_RWS = Mesh(shell, c="red", alpha=0.4)
PSM2_RWS = Mesh(shell, c="green", alpha=0.4)
ECM_FOV = Cone(pos=(0,0,-0.1),r = 0.08 , height = 0.2, res=64, axis=(0,0,1),c="blue",alpha = 0.4)

RCM1xyz = [0.05,0,0]
RCM2xyz = [-0.05,0,0]
PSM1_RWS.pos(RCM1xyz).rotate(30,axis=(0,1,0),point=(RCM1xyz),rad=False)
PSM2_RWS.pos(RCM2xyz).rotate(-30,axis=(0,1,0),point=(RCM2xyz),rad=False)

settings.use_depth_peeling = True
plt = Plotter(shape=(2,1), interactive=False, axes=3)
plt.at(0).show(PSM1_RWS, PSM2_RWS, ECM_FOV, "environment", axes = True)

#intersect
Opt = PSM1_RWS.boolean("intersect", ECM_FOV).c('magenta')
Opt = Opt.boolean("intersect", PSM2_RWS).c('magenta')
print(Opt.volume()*10**6 , "[cc]")
plt.at(1).show(Opt, "intersect volume: %4.2f[cc]" % (Opt.volume()*10**6) , resetcam=False)
plt.interactive().close()
