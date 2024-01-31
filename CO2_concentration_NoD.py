# Use the number of density to calculate the CO2 concentration in the water
# Algorithm:
# 1. Read the number of density file on each timecode
# 2. Normalize the number of density divided by the maximun number of density (necessary?)
# 3. Reorder the number of density from maximun of CO2 number to minimum
# 4. Catch value of CO2 and water from the final frames(ex. final 100 frames) and then average them.
# 5. Calculate the CO2 concentration by the formula: CO2_concentration = CO2_number / (CO2_number + water_number) 


import matplotlib.pyplot as plt
import numpy as np
import ast as ast

# Define the variables
ts = 1000 # timestep of simulation output
dp = 50 # total ts in each frame
tf = 10 # total frames
flag = 'c' # reorder following the certain molecules, 'c' for CO2 and 'w' for water
sf = 0 # starting datas to be averaged (0 is the first frame)
ff = 100 # final datas to be averaged
Concentration = []
frames = np.arange(dp*ts, tf*dp*ts+dp*ts, dp*ts, dtype=int)

for i in range(1,tf+1):
    # 1. Read xvg file
    skip = 24 # skip the first 24 lines
    b = (i-1) * dp * ts  # begining of each frame 
    e = i * dp * ts      # end of each frame
    file_name = "NoD"
    water_file_name = "NoD_water_md_" + str(b) + "_to_" + str(e) + ".xvg"
    CO2_file_name = "NoD_CO2_md_" + str(b) + "_to_" + str(e) + ".xvg"
    data_w = np.genfromtxt(water_file_name, dtype=float, skip_header=skip)
    data_c = np.genfromtxt(CO2_file_name, dtype=float, skip_header=skip)
    slice_number = len(data_w[:,0])

    # 2. Normalize the number of density divided by the maximun number of density (not necessary)
    data_c_normalized = np.zeros((slice_number,2))
    data_w_normalized = np.zeros((slice_number,2))

    max_water = np.max(data_w[:,1])
    max_CO2 = np.max(data_c[:,1])
    max_water_index = np.argmax(data_w[:,1])
    max_CO2_index = np.argmax(data_c[:,1])

    data_w_normalized[:,0] = data_w[:,0]
    data_c_normalized[:,0] = data_c[:,0]
    data_w_normalized[:,1] = data_w[:,1] / max_water
    data_c_normalized[:,1] = data_c[:,1] / max_CO2

    # Create a line plot of the normalized number of density
    plt.plot(data_w[:,0], data_w[:,1], label='water')
    plt.plot(data_c[:,0], data_c[:,1], label='CO2')
    plt.xlabel('Coordinate (nm)')
    plt.ylabel('Number of density (1/nm^3)')
    plt.title('Number of density of water and CO2')
    plt.savefig("Normalized_NoD.png")

    # 3. Reorder the number of density from minimum of CO2 number to maximun
    # Get the indices that would sort the array by the specified column in descending order
    if flag == 'c':
        sorted_indices = np.argsort(data_c[:,1])
    elif flag == 'w':
        sorted_indices = np.argsort(data_w[:,1])
    else:
        print("Please enter the correct flag")
    # Rearrange the array based on the sorted indices
    data_w = data_w[sorted_indices]
    data_c = data_c[sorted_indices]

    # 4. Catch value of CO2 and water from the final frames(ex. final 100 frames) and then average them.
    # 5. Calculate the CO2 concentration by the formula: CO2_concentration = CO2_number / (CO2_number + water_number) 
    Concentration.append( np.mean(data_c[sf:ff,1]) / (np.mean(data_c[sf:ff,1]) + np.mean(data_w[sf:ff,1])))
    print(np.mean(data_c[sf:ff,1]))
    print(np.mean(data_w[sf:ff,1]))
    print(Concentration) # Need to be clarified


# Save the data to a text file
np.savetxt(file_name + ".txt", np.transpose([frames, Concentration]), fmt='%s')
