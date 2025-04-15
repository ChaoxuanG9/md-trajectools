#import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
from matplotlib.pyplot import cm
from time import process_time
from scipy.spatial import distance


files_path = None # pathway to MD trajectory files in xyz format
trj_df = pd.read_csv(files_path, delim_whitespace=True, comment='#')


def pbc2primitive(coordinates, box_size):
    # Apply periodic boundary conditions to handle both positive and negative coordinates
    return (coordinates % box_size + box_size) % box_size

box_size = 14.929668
# Create a new DataFrame with the same columns
trj_prim_df = pd.DataFrame(columns=trj_df.columns)

# Apply periodic boundary conditions to the coordinates using the function
trj_prim_df[['x', 'y', 'z']] = trj_df[['x', 'y', 'z']].apply(lambda col: pbc2primitive(col, box_size))

# Copy the non-coordinate columns to the new DataFrame
trj_prim_df[['timestep', 'element', 'id']] = trj_df[['timestep', 'element', 'id']]

trj_prim_df


def radial_distribut(dataframe, element1, element2, r_max, num_atoms, box_size, num_frames, excluded_frames=0, frame_sample=1, dr=0.2):
    
    bins = int(r_max / dr)
    local_dens_sum = np.zeros(bins, dtype=float)
    local_CN_sum = np.zeros(bins, dtype=float)
    volume = box_size**3

    for frm in range(excluded_frames, num_frames, frame_sample):
        print(frm, end="\r")
        snap = dataframe.iloc[num_atoms * frm: num_atoms * (frm + 1)]

        # Extract positions of the specified elements
        positions_element1 = snap[snap['element'] == element1][['x', 'y', 'z']].to_numpy()

        positions_element2_df = snap[snap['element'] == element2][['x', 'y', 'z']]
        
        distances = np.array([])
        # have PBC included along x and y direction
        index = [0.0, -1.0, 1.0]  # in primitive cell and neighbor cells
        for tran_xx in index:
            for tran_yy in index:
                for tran_zz in index:
                    positions_element2 = positions_element2_df.copy()
                    positions_element2['x'] += tran_xx * box_size
                    positions_element2['y'] += tran_yy * box_size
                    positions_element2['z'] += tran_zz * box_size

                    positions_element2 = positions_element2.to_numpy()

                    # Calculate pairwise distances
                    # The result is a matrix where the element at position (i, j) represents the distance 
                    # between the i-th point in positions_element1 and the j-th point in positions_element2.
                    distances = np.concatenate((distances, distance.cdist(positions_element1, positions_element2, 'euclidean').flatten()))

        # Calculate RDF
        distribution, bin_edges = np.histogram(distances, bins=bins, range=(0, r_max), density=False)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

        
        vol_bins = 4/3 * np.pi * (bin_edges[1:bins+1]**3 - bin_edges[0:bins]**3)
        local_dens = distribution / len(snap[snap['element'] == element1]) / vol_bins
        local_CN = np.cumsum(distribution) / len(snap[snap['element'] == element1])
        
        # local_dens = np.zeros_like(distribution, dtype=float)
        # for bin_count, atom_count in enumerate(distribution):
        #     local_dens[bin_count] = (atom_count / ((4/3)*np.pi*((bin_edges[bin_count+1]**3)-(bin_edges[bin_count]**3))))
             
        local_dens_sum += local_dens
        local_CN_sum += local_CN
    
    rdf = local_dens_sum / len(range(excluded_frames, num_frames, frame_sample)) / (len (snap[snap['element'] == element2]) / volume)
    CN = local_CN_sum / len(range(excluded_frames, num_frames, frame_sample))

    return bin_centers, rdf, CN


bin_centers, rdf, CN = radial_distribut(trj_prim_df, 'N', 'O', 20, 322, 14.929668, 12000, frame_sample=10)



# plt.plot(bin_centers, rdf, label='O-N RDF')
plt.plot(bin_centers, CN / ((4/3) * np.pi * bin_centers**3), label='cumulative local density')
plt.plot(bin_centers, (28 / box_size**3) * np.ones_like(bin_centers), label='bulk density')
plt.xlabel('Distance (Å)')
plt.ylabel('number density (#/Å^3)')
plt.title('N-O')
plt.legend()
plt.show()

