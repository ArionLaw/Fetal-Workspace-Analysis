import vedo
import numpy as np
from ambf.ambf_controller.dvrk.scripts.psmFK import*


PI = np.pi()
# joint max limits: Outer Yaw, Outer Pitch, Insertion, Wrist Roll, Wrist Pitch, Wrist Yaw 
q_max = [90*PI/180 , 54*PI/180 , 0.240 , 180*PI/180 , 180*PI/180 , 180*PI/180 ] # in radians and meters
q_min = [-90*PI/180 , -54*PI/180 , 0.000 , -180*PI/180 , -180*PI/180 , -180*PI/180 ] # in radians and meters

# joint soft limits: Outer Yaw, Outer Pitch, Insertion, Wrist Roll, Wrist Pitch, Wrist Yaw 
qlim_u = [60*PI/180 , 45*PI/180 , 0.18 , 170*PI/180 , 90*PI/180 , 90*PI/180 ] # in radians and meters
qlim_l = [-60*PI/180 , -45*PI/180 , 0.005 , -170*PI/180 , -90*PI/180 , -90*PI/180 ] # in radians and meters

step_r = 2*PI/180 # step size for revolute joints in radians
step_p = 0.005 # step size for prismatic telescopic joint 3 in meters

"""
q1 = np.arange(qlim_l[1],qlim_u[1],step) #Outer Yaw
q2 = np.arange(qlim_l[2],qlim_u[2],step) #Outer Pitch
q3 = np.arange(qlim_l[3],qlim_u[3],step) #Insertion
q4 = np.arange(qlim_l[4],qlim_u[4],step) #Wrist Roll
q5 = np.arange(qlim_l[5],qlim_u[5],step) #Wrist Pitch
q6 = np.arange(qlim_l[6],qlim_u[6],step) #Wrist Yaw
"""

q1 = qlim_l[0]
q2 = qlim_l[1]
q3 = qlim_l[2]
while q1 <= qlim_u[0]:
    while q2 <= qlim_u[1]:
        while q3 <= qlim_u[2]:
            base_frame_pos = np.mat([[0],[0],[0],[0]])
            EE_pos = np.matmul(compute_FK(joint_pos),base_frame_pos)
            print(EE_pos)
            print(compute_FK(joint_pos))
            
            q3 = q3 + step_p
        q2 = q2 + step_r
    q1 = q1 + step_r