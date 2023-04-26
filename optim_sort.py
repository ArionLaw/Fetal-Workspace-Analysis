import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

T = np.load("data_file.npy") #iteration results: output data table
print(T)
np.savetxt("raw-data.csv",T,delimiter=",")
print(T)

# Port to lesion distance
# TIM: is it more important that scope has high dist than tools?
I_score = np.zeros((210,1))
I_best = np.zeros((10,7))
wi1 = 1   #ECM
wi2 = 1   #PSM1
wi3 = 1   #PSM2

# Port to lesion normal angle
# TIM: is it more important that scope be head-on than tools?
J_score = np.zeros((210,1))
J_best = np.zeros((10,7))
wj1 = .2   #ECM
wj2 = .4   #PSM1
wj3 = .4   #PSM2

# Port twist angle (from uterus normal)
# These are all equal as each produces equivalent (undesired) twist
# Different instrument diameters? Maybe smaller tools than scope can get away with more twist?
K_score = np.zeros((210,1))
K_best = np.zeros((10,7))
wk1 = 1/3   #ECM
wk2 = 1/3   #PSM1
wk3 = 1/3   #PSM2

# Scoring each config on each param based on weights
for row in range(210):
     I_score[row,0] = (wi1*T[row,7] + wi2*T[row,8] + wi3*T[row,9]);  # Port - lesion (approach) distance
     J_score[row,0] = (wj1*T[row,10] + wj2*T[row,11] + wj3*T[row,12]);   # Port-axis - lesion-normal (approach) angle
     K_score[row,0] = (wk1*T[row,13] + wk2*T[row,14] + wk3*T[row,15]);     #Port-axis - uterus-normal (twist) angle

# Normalizing scores (clean this up if there's time)
# min value in score
I_min = np.amin(I_score)
J_min = np.amin(J_score)
K_min = np.amin(K_score)

# range of scores
I_range = np.ptp(I_score)
J_range = np.ptp(J_score)
K_range = np.ptp(K_score)

# normalize each param score over 0:1
for row in range(210):
     I_score[row,0] = (I_score[row,0] - I_min) / I_range
     J_score[row, 0] = (J_score[row, 0] - J_min) / J_range
     K_score[row, 0] = (K_score[row, 0] - K_min) / K_range

# Adding normalized scores table to data table
T2 = np.concatenate((T,I_score,J_score,K_score),axis=1)

# Optimizing all 3 scores
overall_ranking = np.zeros((210,1))

wI = .4 # maximize approach distance
wJ = .4 # minimize approach angle
wK = .2 # minimize twist angle

for row in range(210):
     # add max, subtract min => highest overall ranking = optimal solution
     overall_ranking[row,0] = (wI*T2[row,16] - wJ*T2[row,17] - wJ*T2[row,18]);

# #1 optimal ranking + related inputs (fetal tilt and port config)
opt1 = np.max(overall_ranking,0)
opt1_pos = np.argmax(overall_ranking,0)
print('opt1=',opt1,'is at:',opt1_pos)
print('the opt tilt is:',T2[opt1_pos,0])
print('the opt port locs are: ECM(',T2[opt1_pos,1],',',T2[opt1_pos,2],') PSM1(',T2[opt1_pos,3],',',T2[opt1_pos,4],') PSM2(',T2[opt1_pos,5],',',T2[opt1_pos,6],')')

# Finding local maxima (a set of optimized sol'ns)
peaks, _ = find_peaks(overall_ranking[:,0], distance=20)
np.diff(peaks)
set_of_sols = np.zeros((len(peaks),8))
set_of_sols[:,0:6]=(T2[peaks,0:6])
set_of_sols[:,7]=(overall_ranking[peaks,0])
np.savetxt("SetOfSols.csv",set_of_sols,delimiter=",")

# Plot each score (I,J,K on y axis) against each port config (incrementing on x-axis)
xpoints = np.arange(start=0,stop=210,step=1)
plt.xlabel("port config iteration index")
plt.ylabel("optimized score")
plt.plot(xpoints,T2[:,16])
plt.plot(xpoints,T2[:,17])
plt.plot(xpoints,T2[:,18])
plt.plot(xpoints,overall_ranking[:,0],linewidth=4.0)
plt.plot(opt1_pos,opt1,'bo',)
plt.plot(peaks, overall_ranking[peaks], "bx",markersize=7)
plt.title("Port Config Parameter Scoring",loc='left')
plt.legend(["Approach Distance","Approach Angle","Twist Angle","Cost Function","Optimal","Optimal Set"],loc=2,ncol=2,bbox_to_anchor=(0.5,1.3),fancybox=True,shadow=True)
plt.show()

# Sorting based on each "best-of" score
I_sort = T2[T2[:,16].argsort()]
J_sort = T2[T2[:,17].argsort()]
K_sort = T2[T2[:,18].argsort()]

# Producing Top 10 for each param - flipped some params to have best in row 1
I_best = np.delete(np.flipud(I_sort),np.s_[10:211],0)
J_best = np.delete(np.flipud(J_sort),np.s_[10:211],0)
K_best = np.delete(K_sort,np.s_[10:211],0)

# Save results to csv
np.savetxt("MaxLesionDist3.csv",I_best,delimiter=",")
np.savetxt("MinLesionAngle.csv",J_best,delimiter=",")
np.savetxt("MinTwistAngle.csv",K_best,delimiter=",")


#########################################################################################
# QA

# sample array sort to ensure consistent results
     # test_array = np.array([[24,36,3],[1,4,52],[12,32,46],[18,11,23]])
     # print(test_array)
     # test_arrayd = test_array[test_array[:,1].argsort()]
     # print(test_arrayd)

# np.savetxt("foo2.csv", T, delimiter="], [") #output to .csv for visual check
#########################################################################################