import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np

# Load the PDB file
file_name = "md_corrected_test"
pdb_file = file_name + ".pdb"
u = mda.Universe(pdb_file)


# Initialize constants and arrays
nf = u.trajectory.n_frames
natm =  len(u.atoms)
atom_names = u.atoms.names
residue_names = u.atoms.resnames  
guest_index = []
host_index = []
for i in range(0, natm):
    if residue_names[i] == 'CO2' and atom_names[i] == 'C': guest_index.append(i)
    if residue_names[i] == 'SOL' and atom_names[i] == 'OW': host_index.append(i)
nguest = len(guest_index)
nhost = len(host_index)

sf = 0
ef = nf
smooth = 5
yes_mat = np.zeros(nguest)
yes_mat2d = np.zeros((nguest, nguest))
hicnt = np.zeros(nf)
Concentration = np.zeros(nf)

np.savetxt("Cluster.txt", ["Frame  Cluster_number  Concentration"], fmt='%s')

# Iterate through frames
for tf in u.trajectory[sf:ef]: 
    print(tf)
    cf = tf.frame # Current frame
    for ts in u.trajectory[tf.frame: tf.frame+smooth]:
        print(ts)
        # Access information for each frame
        coordinates = u.atoms.positions
        atom_names = u.atoms.names
        residue_names = u.atoms.resnames  
    
        # Record the molecular index in the cutoff range
        # The first pick in RDF of liquid CO2 is 4.0 A.
        # Here takes 4.2 A as cut-off distance
        for i in range(0, nguest):
            for j in range(i+1, nguest):
                # Calculate the distance between two carbon atoms
                if np.linalg.norm(coordinates[guest_index[i]]-coordinates[guest_index[j]]) < 4.2:
                    yes_mat2d[i, j] += 1
                    yes_mat2d[j, i] += 1                    
    
    # After ts loop, ts would become 0.(Bug?) So ts have to be turned back.
    ts = u.trajectory[cf] 
    for i in range(0, nguest):
        for j in range(i + 1, nguest):
            if yes_mat2d[i, j] == smooth:
                yes_mat[i] = 1
                yes_mat[j] = 1
    
    icount = 0
    for i in range(nguest):
        if yes_mat[i] == 1:
            icount += 1

    hicnt[cf] = icount
    Concentration[cf] = (nguest-hicnt[cf])/(nguest+nhost)
    yes_mat2d[:, :] = 0
    yes_mat[:] = 0
    print("Current frame={:d}, Cluster_number={:d}, Concentration={:+.4f} \n".format(cf, icount, Concentration[cf]))
    
    # Output to the file
    record = "{:d}  {:d}  {:+.4f} \n".format(cf, icount, Concentration[cf])
    with open("Cluster.txt", 'a') as f:
        f.write(record)

# Open the file in read mode
with open("Cluster.txt", 'r') as file:
    # Read all lines from the file
    lines = file.readlines()[1:]

# Initialize empty arrays to store data
time = []
Xco2 = []

# Skip the first line (header)
lines = lines[1:]

# Iterate over each line and extract values
for line in lines:
    # Split the line into columns
    columns = line.split()

    # Append values to the arrays
    time.append(float(columns[0])*200)
    Xco2.append(float(columns[2]))

# Create a line plot
plt.plot(time, Xco2)

# Add labels and a title
# plt.xlabel('Time (ns)')
plt.xlabel('Time(ps)')
plt.ylabel('Concentration(mole fraction)')
plt.title('CO2 Concentration-time plot')

# Display the plot (optional)
# plt.show()
plt.savefig("Cluster.png")