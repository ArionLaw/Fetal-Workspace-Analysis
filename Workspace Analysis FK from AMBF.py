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
qlim_u = [30*PI/180 , 30*PI/180 , 0.18 , 120*PI/180 , 90*PI/180 , 90*PI/180 ] # in radians and meters
qlim_l = [-30*PI/180 , -30*PI/180 , 0.04 , -120*PI/180 , -90*PI/180 , -90*PI/180 ] # in radians and meters

step_r_L1 = 5*PI/180 # coarse step size for revolute joints in radians (used for mesh refinement)
step_r_L2 = 2.5*PI/180 # coarse step size for revolute joints in radians (used for mesh refinement)
step_r_L3 = 2*PI/180 # fine step size for revolute joints in radians (used for mesh refinement)
step_r_L4 = 1*PI/180 # fine step size for revolute joints in radians (used for mesh refinement)
step_p = 0.005 # step size for prismatic telescopic joint 3 in meters

q1 = qlim_l[0]
q2 = qlim_l[1]
q3 = qlim_l[2]

i = 0
PSM_EE_xyz = [[]]
tol = 0.006 #numerical operation tolerance
while q3 <= qlim_u[2]:
    while q2 <= qlim_u[1]+tol:
        while q1 <= qlim_u[1]+tol:
            q1 = round(q1,3)
            q2 = round(q2,3)
            q3 = round(q3,3)

            joint_pos = [q1 , q2 , q3 , 0 , 0 , 0]
            base_frame_pos = np.mat([[0],[0],[0],[1]])
            PSM_EE_pos = np.matmul(compute_FK(joint_pos),base_frame_pos)
            #print("Outer Yaw [rads]: " , q1)
            #print("Outer Pitch [rads]: " , q2)
            #print("Insertion [m]: " , q3 )
            #print("End Effector Position[m]: ")
            #print( PSM_EE_pos , '\n')"""
            PSM_EE_xyz.append([float(PSM_EE_pos[0]),float(PSM_EE_pos[1]),float(PSM_EE_pos[2])])
            i = i + 1

            if q3 <= qlim_l[2] + step_p:
                q1 = q1 + step_r_L1
            if q3 <=  qlim_u[2] - 20*step_p:
                q1 = q1 + step_r_L2
            if q3 <= qlim_u[2] - 1*step_p:
                q1 = q1 + step_r_L3
            else: 
                q1 = q1 + step_r_L4
        
        q1 = qlim_l[0]
        if q3 <= qlim_l[2] + step_p:
            q2 = q2 + step_r_L1
        if q3 <= qlim_u[2] - 20*step_p:
            q2 = q2 + step_r_L2
        if q3 <= qlim_u[2] - 1*step_p:
            q2 = q2 + step_r_L3
        else:
            q2 = q2 + step_r_L4
    
    print("Insertion [m]: " , q3 )
    q2 = qlim_l[1]
    q3 = q3 + step_p

PSM_EE_xyz.pop(0)
#print(PSM_EE_xyz)

cloud = pv.PolyData(PSM_EE_xyz)
cloud.plot()
volume = cloud.delaunay_3d(alpha = 0.004)
shell = volume.extract_geometry()
#shell.plot()

PSM1_RWS = Mesh(shell, c="red", alpha=0.6)
show(PSM1_RWS, axes = True)   