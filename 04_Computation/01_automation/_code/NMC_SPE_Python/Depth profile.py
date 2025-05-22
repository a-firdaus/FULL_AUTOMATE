import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.ticker import MaxNLocator
import os
import tkinter as tk
from tkinter import filedialog
from vpython import vector as vpythvector
from molatomclasses import atom

def add_atomi_to_molecule(molecule, allatompos):
    for atomid, atomposlocation in allatompos.items():
        molecule.append(atom(atomid, ''.join([i for i in atomid if not i.isdigit()]), atomposlocation))
def create_atomsymbolfulllist(atoms, quantities):
    atomsymbolfulllist = []
    for index_ in range(len(atoms)):
        for j in range(quantities[index_]):
            atomsymbolfulllist.append(atoms[index_])
    return atomsymbolfulllist
def select_file():
    """Opens a file dialog to select a file and returns the file path."""
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_path = filedialog.askopenfilename(title="Select a file")
    if not file_path:
        raise FileNotFoundError("No file selected.")
    return file_path
def read_file(file_path):
    """Reads a file and returns its content as a list of lines."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines
def treat_lines(lines):
    """Processes lines to extract element symbols and positions.

    Args:
        lines (list): List of strings, each representing a line from the file.

    Returns:
        tuple: A list of element symbols and a numpy array of positions.
    """

    atoms = []
    quantities = []
    history = []
    nbr_line = -9
    check = False
    avector = None
    bvector = None
    cvector = None
    latticevectorarray = None
    currentindex = -1
    allatompos = {}
    molecule = []

    linecount = 0
    for line in lines:
        nbr_line += 1
        linecount += 1
        if line.strip() == "Selective dynamics":
            check = True
            atoms = history[-2].split()
            quantities = [int(x) for x in history[-1].split()]
            atomsymbolfulllist = create_atomsymbolfulllist(atoms, quantities)

        if not check:
            history.append(line.strip())

        if linecount == 3:
            avector = np.array([float(k) for k in line.strip().split()])
            continue
        elif linecount == 4:
            bvector = np.array([float(k) for k in line.strip().split()])
            continue
        elif linecount == 5:
            cvector = np.array([float(k) for k in line.strip().split()])
            latticevectorarray = np.stack([avector, bvector, cvector], axis=0)
            continue
        elif linecount > 9 and linecount < sum(quantities, 10):
            currentindex += 1
            pos = line.strip().split()[:3]
            #location_ = vpythvector(*[float(pos[0]), float(pos[1]), float(pos[2])])
            location_ = [float(pos[0]) * np.linalg.norm(avector), float(pos[1]) * np.linalg.norm(bvector), float(pos[2]) * np.linalg.norm(cvector)]
            currentsymbol = atomsymbolfulllist[currentindex] + str(currentindex)
            allatompos[currentsymbol] = location_

    add_atomi_to_molecule(molecule, allatompos)

    return molecule, latticevectorarray, quantities, atomsymbolfulllist
def exclude_spe(spe, atomsymbolfulllist, comp_spe):
    spe_elem = ['C', 'H', 'O', 'N']
    electrode_selection_atoms = []
    electrode_selection_ids = []

    for i, atom in enumerate(atomsymbolfulllist):
        if atom not in spe_elem:
            electrode_selection_ids.append(i)
        else:
            if atom == 'O' and spe in ['Alcohol', 'Ester', 'Carbonate', 'Ether', 'Urethane']:
                electrode_selection_ids.append(i)

    numberO = comp_spe[spe]['O']
    oxygen_indices = [i for i in electrode_selection_ids if atomsymbolfulllist[i] == 'O']
    excluded_oxygen_indices = set(oxygen_indices[-numberO:])
    final_indices = [i for i in electrode_selection_ids if i not in excluded_oxygen_indices]
    electrode_selection_atoms = [atomsymbolfulllist[i] for i in final_indices]

    return electrode_selection_atoms, final_indices
def plot_depth_profile(atoms):

    all_positions = np.array([atom.pos for atom in atoms])
    z_positions = all_positions[:, 2]  # Extract z-coordinates
    atomsymbolfulllist = [atom.symbol for atom in atoms]
    unique_elements = sorted(set(atomsymbolfulllist))
    z_lim = max(z_positions)

    # Define bins for z-axis
    bins = np.linspace(0, round(max(z_positions),2)+0.01, num=10) # Adjust 'num=' as needed
    cumulative_composition = {element: [] for element in unique_elements}

    for i in range(len(bins) - 1):
        # Find indices of atoms within the current bin
        in_bin = (z_positions >= bins[i]) & (z_positions < bins[i + 1])
        total_in_bin = np.sum(in_bin)
        print('-'*60)
        print(f"There are {total_in_bin} atoms in this bin ({bins[i]} - {bins[i+1]})")
        for element in unique_elements:
            element_count = sum(e == element for e, inside in zip(atomsymbolfulllist, in_bin) if inside)
            print(f"There are {element_count} {element} atoms in this bin")
            cumulative_composition[element].append(element_count / len(atomsymbolfulllist) * 100)

            # if total_in_bin > 0 or i == 0:
            #      cumulative_composition[element].append(element_count/len(atomsymbolfulllist)*100)
            # else:
            #      cumulative_composition[element].append(cumulative_composition[element][-1])


    # Plot the depth profile
    fig, ax = plt.subplots(figsize=(12,7))
    for element in unique_elements:
        ax.plot(bins[:-1], cumulative_composition[element], label=element)#, marker='o')
        #ax.bar(bins[:-1], cumulative_composition[element], label=element, width=np.diff(bins)[0], align='edge', alpha=0.7)
    ax.set_xlabel("Depth (z-axis)")
    ax.set_ylabel("Atomistic composition (% of atoms)")
    ax.set_xlim(0, z_lim)
    ax.set_ylim(0, 20)
    ax.legend()
    ax.grid(False)
    plt.title("Depth Profile Composition")
    plt.show()

def main():
    from molatomclasses import atom

    comp_spe = {'Alcohol':{'O':3},
                'Amine':{'N':3},
                'Carbonate':{'O':6},
                'Ester':{'O':4},
                'Ether':{'O':2},
                'Urethane':{'O':4, 'N':2}} # Check if this is correct if any issues
    try:
        spe = input("What is the SPE on top of the electrode?\n")
        while spe not in comp_spe.keys():
            print("This is not a valid input!")
            spe = input("What is the SPE on top of the electrode?\n")
        file_path = select_file()
        lines = read_file(file_path)
        atoms, latticevectorarray, quantities, atomsymbolfulllist = treat_lines(lines)

        electrode_selection_atoms, electrode_selection_ids = exclude_spe(spe, atomsymbolfulllist, comp_spe)
        electrode_only_atoms = []
        for at, id in zip(electrode_selection_atoms, electrode_selection_ids):
            electrode_only_atoms.append(atom(id, at, atoms[id].pos))

        elem_counts_dict = {}
        for i, elem in enumerate(sorted(set(atomsymbolfulllist))):
            elem_counts_dict[elem] = quantities[i]
        plot_depth_profile(electrode_only_atoms)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()