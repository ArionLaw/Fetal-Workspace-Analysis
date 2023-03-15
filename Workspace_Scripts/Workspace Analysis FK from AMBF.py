import vedo
import numpy as np
import sys

from Fetal_Workspace_Analysis.ambf.ambf_controller.dvrk.scripts.psmFK import*

# joint limits
q1_lim = [0 , 90]
q2_lim = [0 , 90]
q3_lim = [0 , 90]
q4_lim = [0 , 90]
q5_lim = [0 , 90]
q6_lim = [0 , 90]
step = 5 #rad

joint_pos = [0 , 0 , 0 , 0 , 0 , 0]
base_frame_pos = np.mat([[0],[0],[0],[0]])
EE_pos = np.matmul(compute_FK(joint_pos),base_frame_pos)
print(EE_pos)
print(compute_FK(joint_pos))

