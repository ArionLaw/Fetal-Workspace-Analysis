import vedo
import numpy as np
import kinpy as kp
import math

chain = kp.build_chain_from_urdf(open("From Github/kinpy/examples/kuka_iiwa/model.urdf").read())
print(chain)
# lbr_iiwa_link_0_frame
# └──── lbr_iiwa_link_1_frame
#       └──── lbr_iiwa_link_2_frame
#             └──── lbr_iiwa_link_3_frame
#                   └──── lbr_iiwa_link_4_frame
#                         └──── lbr_iiwa_link_5_frame
#                               └──── lbr_iiwa_link_6_frame
#                                     └──── lbr_iiwa_link_7_frame

print(chain.get_joint_parameter_names())
# ['lbr_iiwa_joint_1', 'lbr_iiwa_joint_2', 'lbr_iiwa_joint_3', 'lbr_iiwa_joint_4', 'lbr_iiwa_joint_5', 'lbr_iiwa_joint_6', 'lbr_iiwa_joint_7']

th = {'lbr_iiwa_joint_2': math.pi / 4.0, 'lbr_iiwa_joint_4': math.pi / 2.0}
ret = chain.forward_kinematics(th)