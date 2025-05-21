from ase import Atoms
import numpy as np

from pymatgen.core import Molecule

def align_principal_axes(atoms: Atoms):
    pos = atoms.get_positions()
    pos -= pos.mean(axis=0)  # Center molecule

    # Compute inertia tensor
    inertia = np.dot(pos.T, pos)
    eigvals, eigvecs = np.linalg.eigh(inertia)

    # Align molecule along principal axes
    new_pos = np.dot(pos, eigvecs)
    atoms.set_positions(new_pos)
    return atoms

def substitute_ni_atoms(structure, nmc_ratio, seed=None):
    nmc_structure = structure.copy()
    
    nmc_ratio = np.array(nmc_ratio)
    nmc_ratio = nmc_ratio/nmc_ratio.sum()
    
    # Use a NumPy RNG like QuantumATK does
    rng = np.random.RandomState(seed)

    # Find all Ni atoms
    ni_indices = [i for i, site in enumerate(nmc_structure) if site.specie.symbol == "Ni"]
    print("Original Ni indices:", ni_indices)
    
    n_total = len(ni_indices)
    
    # Compute target numbers for each element
    ni_frac, mn_frac, co_frac = nmc_ratio
    n_ni = round(n_total * ni_frac)
    n_mn = round(n_total * mn_frac)
    n_co = n_total - n_ni - n_mn  # Ensure total is preserved
    
    # Shuffle Ni indices
    rng.shuffle(ni_indices)
    print("Shuffled Ni indices:", ni_indices)
    
    # Replace Ni with Mn
    for i in ni_indices[n_ni:n_ni + n_mn]:
        nmc_structure.replace(i, "Mn")
    
    # Replace Ni with Co
    for i in ni_indices[n_ni + n_mn:n_ni + n_mn + n_co]:
        nmc_structure.replace(i, "Co")
    
    print(f"NMC ratio: {nmc_ratio}")
    print(f"n_ni: {n_ni}, n_mn: {n_mn}, n_co: {n_co}")

    nmc_structure = nmc_structure.get_sorted_structure()
    
    return nmc_structure

def get_principal_axis(mol: Molecule):
    coords = mol.cart_coords - mol.center_of_mass
    cov = np.cov(coords.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # Return the eigenvector with the largest eigenvalue (longest axis)
    return eigvecs[:, np.argmax(eigvals)]