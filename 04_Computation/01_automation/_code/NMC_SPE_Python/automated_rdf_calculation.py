import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from tkinter import Tk, filedialog, simpledialog
from typing import Dict, Optional
from molatomclasses import atom
from vpython import vector as vpythvector
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib.lines import Line2D
from pyautogui import confirm
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.topology.guessers import guess_bonds, guess_angles
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.core.groups import AtomGroup
from scipy.ndimage import gaussian_filter1d
from units import unit
matplotlib.use('TkAgg')


radiuslist={"H":1.20,"Li":1.82,"Cl":2.75 ,"C":1.70, "Si":2.10 ,"N":1.55, "O":1.52,"P":1.80,"S":1.80,'Ne':1.54,'Ar':1.88,
            'He':1.40,"Na":2.27,"F":1.47,"Br":1.85,"Tc":2.05,"Pm":2.43,"Cd":1.55,"Co":1.35,"Mn":1.61,"Fe":1.61,"Ni":1.49}
colourdict = {"H": (0.8, 0.8, 0.8), "Li":(252/255, 144/255, 245/255),"Cl": (0, 1, 0), "C": (0, 0, 0), "O": (1, 0, 0), "N": (0, 0, 1), "Ar": (0.90, 0.90, 0.90),
              "Ne": (0.90, 0.90, 0.90), "S": (1, 1, 0), "He": (0.90, 0.90, 0.90), "Na": (1, 1, 1), "Pm": (1, 0.8, 1),
              "Tc": (0.8, 1, 1), "F": (159 / 255.0, 252 / 255.0, 206 / 255.0), "Br": (153 / 255.0, 34 / 255.0, 0),
              "Cd": (0.8, 1, 1), "Co": (0.8, 1, 1),"Si":(153/255, 51/255, 242/255),"Ni":(0.5, 1, 1),"Mn":(0.3, 1, 1),"P":(180/255.0, 52/255.0, 235/255.0),"Fe":(0.3, 1, 1)}
molarweightlist={"H":1.00794,"O":15.9994,"C":12.0107,"Ar":39.948,"He":4.002602,"Ne":20.1797,"N":14.0067,"Cl":35.453,
                 "Na":22.989769,"F":18.998403,"Br":79.904,"Tc":98.0,"Pm":145.0,"Cd":112.411,"Co":58.9332,"Pb":207.2,
                 'Cu':63.546,"Fe":55.845}
element_colors = {'H': 'oldlace', 'Li': 'green', 'O': 'red', 'C': 'saddlebrown', 'N': 'blue',
                      'Mn': 'darkorchid', 'Ni': 'silver', 'Co': 'indigo'}
comp_spe = {'Alcohol':{'O':3},
            'Amine':{'N':3},
            'Carbonate':{'O':6},
            'Ester':{'O':6},
            'Ether':{'O':3},
            'Urethane':{'O':4, 'N':2}}


def select_folder() -> Optional[str]:
    """
    Opens a Tkinter dialog for the user to manually select a folder containing XDATCAR files.

    Returns:
        The path to the selected folder, or None if the selection fails.
    """
    # Initialize Tkinter root and hide the main window
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)  # Bring dialog to the front

    folder_path = filedialog.askdirectory(title="Select Folder Containing XDATCAR Files")
    root.destroy()  # Close the Tkinter root window

    return folder_path if os.path.isdir(folder_path) else None
def get_xdatcar_files(folder_path: str) -> Dict[str, str]:
    """
    Retrieves all files in the specified folder that start with 'XDATCAR_'.

    Args:
        folder_path (str): Path to the directory containing the files.

    Returns:
        dict: A dictionary with file names as keys and their full paths as values.
    """
    xdatcar_files = {
        filename: os.path.join(folder_path, filename)
        for filename in os.listdir(folder_path)
        if filename.startswith("XDATCAR_")
    }
    return xdatcar_files
def parse_xdatcar_name(file_name: str) -> Dict[str, str]:
    """
    Parses the XDATCAR file name to build a hierarchical structure.

    Args:
        file_name (str): The name of the XDATCAR file (e.g., 'XDATCAR_LiNMC622_Alcohol').

    Returns:
        dict: A nested dictionary structure with keys based on the parsed name.
    """
    parts = file_name.split('_')

    if len(parts) < 2:
        return {}  # Skip files that don't fit the expected format

    _, li_type, material, solvent = parts[0], parts[1][:2], parts[1][3:], parts[2]

    structure = {
        'With Li' if li_type == 'Li' else 'Without Li': {
            material: {solvent: ""}
        }
    }
    return structure
def add_paths_to_dict(xdatcar_files: Dict[str, str]) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Builds the hierarchical dictionary and adds paths to each XDATCAR entry.

    Args:
        xdatcar_files (dict): Dictionary of XDATCAR files and their paths.

    Returns:
        dict: A nested dictionary with paths added based on file name parsing.
    """
    xdatcar_dict = {'With Li': {}, 'Without Li': {}}

    for file_name, file_path in xdatcar_files.items():
        structure = parse_xdatcar_name(file_name)

        if not structure:
            continue

        for li_key, materials in structure.items():
            for material, solvents in materials.items():
                for solvent, _ in solvents.items():
                    # Assign path to the final dictionary level
                    xdatcar_dict.setdefault(li_key, {}).setdefault(material, {})[solvent] = file_path

    return xdatcar_dict
def read_xdatcar(file_path):
    """
    Reads an XDATCAR file and returns its lines as a list.

    Args:
        file_path (str): Path to the XDATCAR file.

    Returns:
        list: List of lines from the XDATCAR file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines
def write_lattice(name, lattice):
    """
    Writes lattice vectors to a file with the given name.

    Args:
        name (str): Name of the output file (without extension).
        lattice (list of list): List of lattice vectors.

    Writes:
        A file with each lattice vector on a new line.
    """
    filename = name + ".vectors"
    with open(filename, "w") as file:
        for sublist in lattice:
            line = " ".join(str(element) for element in sublist)
            file.write(line + "\n")

def system(lines):
    """
    Parses lattice and atom information from XDATCAR lines.

    Args:
        lines (list): List of lines from the XDATCAR file.

    Returns:
        tuple: Contains lattice lengths (a, b, c), full atom symbol list,
               atom quantities, and lattice vectors.
    """

    avector = None
    bvector = None
    cvector = None
    atoms = []
    quantities = []
    lattice = []

    for i, line in enumerate(lines):
        if i == 2:
            lattice.append(line.split())
            avector = np.array([float(k) for k in line.strip().split()]) / 10
            continue
        elif i == 3:
            lattice.append(line.split())
            bvector = np.array([float(k) for k in line.strip().split()]) / 10
            continue
        elif i == 4:
            lattice.append(line.split())
            cvector = np.array([float(k) for k in line.strip().split()]) / 10
            continue
        elif i == 5:
            atoms = line.strip().split()
            continue
        elif i == 6:
            quantities = [int(x) for x in line.strip().split()]
            atomsymbolfulllist = create_atomsymbolfulllist(atoms, quantities)
            break

    a = np.linalg.norm(avector) * 10
    b = np.linalg.norm(bvector) * 10
    c = np.linalg.norm(cvector) * 10
    volume = a * b * c
    print("a length", a, "")
    print("b length", b, "")
    print("c length", c, "")

    return a, b, c, atomsymbolfulllist, quantities, lattice
def create_atomsymbolfulllist(atoms, quantities):
    """
    Creates a list of atom symbols based on the quantities for each element.

    Args:
        atoms (list of str): List of atom types.
        quantities (list of int): List of quantities for each atom type.

    Returns:
        list: Full list of atom symbols.
    """
    atomsymbolfulllist = []
    for index_ in range(len(atoms)):
        for j in range(quantities[index_]):
            atomsymbolfulllist.append(atoms[index_])
    return atomsymbolfulllist
def write_xyz(lines, output_file, a, b, c, atomsymbolfulllist, quantities):
    """
        Writes atom coordinates from XDATCAR lines to an XYZ file and generates an interface of atom objects.

        Args:
            lines (list): List of lines from the XDATCAR file.
            output_file (str): Name of the output XYZ file.
            a (float): Lattice parameter along the x-axis.
            b (float): Lattice parameter along the y-axis.
            c (float): Lattice parameter along the z-axis.
            atomsymbolfulllist (list of str): Full list of atom symbols for each atom in the system.
            quantities (list of int): List containing the quantity of each atom type.

        Returns:
            list: A list of `atom` objects representing the interface, each initialized with an ID, symbol, and location vector.

        Writes:
            A formatted XYZ file with atom types and coordinates scaled by lattice parameters.

        Note:
            In the first configuration (defined by the XDATCAR file), this function constructs `atom` objects with unique IDs
            and coordinates, storing them in the `interface` list for further use.
    """

    with open(output_file, "w") as file:
        first_config = True
        interface = []

        counter = 0
        for i, line in enumerate(lines):
            if "Direct configuration=" in line:
                file.write(str(sum(quantities)) + '\n\n')
                counter = 0
                if i > 7:
                    first_config = False
            elif i > 7:
                x = str('{:.6f}'.format(round(float(line.strip().split()[0]) * a, 6)))
                y = str('{:.6f}'.format(round(float(line.strip().split()[1]) * b, 6)))
                z = str('{:.6f}'.format(round(float(line.strip().split()[2]) * c, 6)))

                elem = atomsymbolfulllist[counter]

                if first_config:
                    location_ = vpythvector(*[float(x), float(y), float(z)])
                    id = str(counter)
                    atomid = elem + id
                    interface.append(atom(atomid, elem, location_))

                counter += 1
                file.write(f"{' ' * (2 - len(elem)) + elem}{' ' * (11 - len(x))}{x + ' ' * (11 - len(y))}{y + ' ' * (11 - len(z))}{z}\n")

    return interface

def draw(interface, uniq_elem, radiuslist, element_colors):
    """
        Draws a 2D scatter plot of atomic positions with customizable size and color
        based on element type and user-selected orientation.

        This function prompts the user to select an axis orientation for the plot
        (XY, YZ, or XZ plane) and uses this choice to set up the 2D view. Each atom
        is represented by a point with size and color corresponding to its element,
        using values from the radius list and element colors dictionary. A legend
        for element types is also added.

        Args:
            interface (list): List of atom objects, each with positional attributes `pos.x`, `pos.y`, and `pos.z`.
            uniq_elem (list of str): Unique elements present in the interface, used for the plot legend.
            radiuslist (dict): Dictionary mapping element symbols to atomic radii for determining point sizes.
            element_colors (dict): Dictionary mapping element symbols to color codes for the scatter plot.

        Returns:
            tuple: Returns a tuple containing:
                - `ax` (matplotlib.axes.Axes): Axes object for the scatter plot.
                - `fig` (matplotlib.figure.Figure): Figure object containing the plot.
                - `pts` (matplotlib.collections.PathCollection): Collection of scatter plot points for interaction.

        Raises:
            NotImplementedError: If an unrecognized orientation is selected.

        Note:
            The user is prompted to select an orientation for the scatter plot. If the orientation is incorrect,
            they can cancel to reselect. After orientation selection, the user can select molecules by clicking,
            and then press enter to finalize the selection or close the window to abort.
        """

    orientation_option = confirm(
        text="Select the appropriate axis view.\nIf this is not the right choice, then press cancel in a later dialog to return here.",
        title="Select the axis orientation", buttons=["XY", "YZ", "XZ"])
    data_firstaxis = []
    data_secondaxis = []
    sizes = []
    colours = []
    for atomi in interface:
        if orientation_option == "XY":
            data_firstaxis.append(atomi.pos.x)
            data_secondaxis.append(atomi.pos.y)
            sizes.append(radiuslist[atomi.symbol] * 20)
            colours.append(element_colors[atomi.symbol])
        elif orientation_option == "YZ":
            data_firstaxis.append(atomi.pos.z)
            data_secondaxis.append(atomi.pos.y)
            sizes.append(radiuslist[atomi.symbol] * 20)
            colours.append(element_colors[atomi.symbol])
        elif orientation_option == "XZ":
            data_firstaxis.append(atomi.pos.x)
            data_secondaxis.append(atomi.pos.z)
            sizes.append(radiuslist[atomi.symbol] * 20)
            colours.append(element_colors[atomi.symbol])
        else:
            raise NotImplementedError
    ax = plt.gca()
    fig = plt.gcf()
    pts = ax.scatter(data_firstaxis, data_secondaxis, color=colours, s=sizes)

    # Lists to store legend handles and labels
    marker_handles = []
    marker_labels = []

    # Add legend handles and labels for markers
    for element in uniq_elem:
        if element not in marker_labels:
            marker_handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=element_colors[element],
                                         markersize=8, label=element))  # Adjust the markersize value as needed
            marker_labels.append(element)
        # Create custom legends for materials and markers
        marker_legend = plt.legend(handles=marker_handles, labels=marker_labels, title='Element', loc='upper left',
                                   bbox_to_anchor=(1.0, 1), fontsize=10, markerscale=1.5, title_fontsize=12)
    ax.set_title("Select the molecule and press enter to accept selected points. Otherwise, close the window.")
    return ax, fig, pts
def accept(event, ax, fig, selector, all_selections):
    if event.key == "enter":
        print("Point selection OK!\nHere are the indices:")
        #all_selections.append(selector.ind.tolist()[0])
        all_selections.extend(selector.ind)
        print(all_selections)
        selector.disconnect()
        ax.set_title("(Selection OK, please close the window)")
        fig.canvas.draw()
def manage_selection_plot(selection_type, interface, uniq_elem, radiuslist, element_colors):
    """
        Manages the selection of atoms in a plot by allowing the user to interactively select points and confirm the selection.

        Args:
            selection_type (str): Type of selection, typically 'Electrode' or another category, which affects confirmation prompts.
            interface (list): List of `atom` objects with positions and element symbols.
            uniq_elem (list): List of unique element symbols present in the interface.
            radiuslist (dict): Dictionary mapping element symbols to atomic radii for point sizing in the plot.
            element_colors (dict): Dictionary mapping element symbols to color codes for plot coloring.

        Returns:
            list: List of selected points after user confirmation, representing atoms chosen from the plot.

        User Interaction:
            - For `Electrode` type, the user confirms selection with an option to retry.
            - For other types, the user is prompted to continue selection or finish based on their choices.

        Process:
            1. Creates an interactive scatter plot by calling `draw()`, displaying atoms with sizes and colors according to element.
            2. A `SelectFromCollection` object is used for atom selection.
            3. User confirmation dialogs allow for iteration until selections are finalized or confirmed.

        Note:
            - Each iteration shows the plot, and selected atoms are stored in `all_selections`.
            - The process repeats until the user confirms a satisfactory selection.
        """

    has_been_selected = False
    all_selections = []

    if selection_type == 'Electrode':
        while not has_been_selected:
            ax, fig, pts = draw(interface, uniq_elem, radiuslist, element_colors)
            fig.canvas.get_tk_widget().focus_set()
            selector = SelectFromCollection(ax, pts)
            fig.canvas.mpl_connect("key_press_event", lambda event: accept(event, ax, fig, selector, all_selections))
            ax.set_title("Select electrode outer layer and press enter.\nOtherwise, close the window.")
            plt.show()
            plt.close()

            if confirm(text="Did you manage to select the atoms?", buttons=["Yes", "No"]) == "No":
                all_selections = []
                continue
            else:
                return all_selections

    else:
        while not has_been_selected:
            ax, fig, pts = draw(interface, uniq_elem, radiuslist, element_colors)
            fig.canvas.get_tk_widget().focus_set()
            selector = SelectFromCollection(ax, pts)
            fig.canvas.mpl_connect("key_press_event", lambda event: accept(event, ax, fig, selector, all_selections))
            ax.set_title("Select entire SPE and press enter.\nOtherwise, close the window.")
            plt.show()
            plt.close()

            if confirm(text="Do you still need to select atoms?", buttons=["Yes", "No"]) == "Yes":
                continue
            else:
            #     if len(all_selections) != sum(comp_spe.get(selection_type, {}).values()):
            #         print("Something went wrong with the selection. Please try again from scratch.")
            #         all_selections = []
            #         continue

                return all_selections
def rdf_calculation(xyz_path, x, y, z, ind_ref, ind, start=1250, nbins = 75, range = (0.0, 6.0)):
    """
    Calculates the radial distribution function (RDF) for a specified pair of atoms in an XYZ file.

    Args:
        xyz_path (str): Path to the XYZ file containing the atomic structure.
        x (float): Length of the simulation box in the x-axis (in Ångströms).
        y (float): Length of the simulation box in the y-axis (in Ångströms).
        z (float): Length of the simulation box in the z-axis (in Ångströms).
        ind_ref (int): Index of the reference atom in the Universe (0-based indexing).
        ind (int): Index of the atom for which to calculate the RDF with respect to the reference atom (0-based indexing).

    Returns:
        InterRDF: MDAnalysis InterRDF object containing RDF data.
                  The RDF can be accessed with `rdf.rdf`, and the bins (distances) with `rdf.bins`.

    """

    u = mda.Universe(xyz_path, format="XYZ")
    u.dimensions = [x, y, z, 90, 90, 90]

    atom1 = AtomGroup([ind_ref], u)
    atom2 = AtomGroup([ind], u)

    rdf = InterRDF(atom1, atom2,
                   nbins = nbins,
                   range = range,
                   )
    rdf.run(start=start)

    return rdf
def rdf_meets_criterion(rdf_obj, criterion_range=4):
    """
    Checks if an RDF (radial distribution function) meets the criterion of having
    zero values for all RDF points within a specified range of distances.

    Args:
        rdf_obj (MDAnalysis.analysis.rdf.InterRDF): An RDF object containing
            computed RDF data and bins.
        criterion_range (float, optional): The maximum distance (in Å) within
            which the RDF values are checked to be zero. Default is 5 Å.

    Returns:
        bool: True if all RDF values within the criterion range are zero, False otherwise.

    Example:
        rdf = InterRDF(atom_group1, atom_group2, nbins=75, range=(0.0, 10.0))
        rdf.run()
        meets_criterion = rdf_meets_criterion(rdf, criterion_range=5)
    """

    # Select bins within the criterion range
    bins_within_range = rdf_obj.bins[rdf_obj.bins <= criterion_range]
    rdf_values_within_range = rdf_obj.rdf[:len(bins_within_range)]

    # Check if all rdf values within the range are 0
    return np.all(rdf_values_within_range == 0)
def plot_rdfs(data, save_path=None, interface_name=None):
    """
    Plots the RDF (radial distribution function) data from a given dictionary
    as subplots with a maximum of three subplots per row. Titles and legends
    are displayed for each subplot, while axis labels are customized based on
    their positions. Optionally saves the generated plot to a specified path.

    Args:
        data (dict): Dictionary containing RDF data for plotting. Expected format:
            {'atom_id': {'label': rdf_object}}, where `atom_id` is a string
            identifier for the atom, `label` is a string representing the specific
            atom pair, and `rdf_object` contains the RDF results.
        save_path (str, optional): File path to save the plot. If provided, the
            plot will be saved to this path instead of being displayed.
        interface_name (str, optional): Combination of electrode and spe material.
            If provided, this name will be used to save the plot instead of displaying it.

    Example of `data`:
        {'C0': {'C0-Li25': rdf}, 'C1': {'C1-Li26': rdf}}

    Notes:
        - Each subplot has its own title (customizable within the function).
        - Legends are visible for each subplot.
        - X-axis is labeled as 'Distance (Å)' with range 0-6 Å, but labels and ticks
          are only displayed on the bottom row if multiple rows are present.
        - Y-axis is labeled as 'g(r)' only on the leftmost subplots of each row.
        - No grid is applied to any subplot.

    Returns:
        None: The function generates and displays the subplots or saves them to `save_path` with
        'interface_name' as file name.
    """

    # Calculate number of subplots and rows
    num_plots = len(data)
    cols = 3
    rows = (num_plots + cols - 1) // cols

    # Set up the subplot grid
    fig, axs = plt.subplots(rows, cols, figsize=(15, 5 * rows), constrained_layout=True)
    axs = axs.flatten()  # Flatten for easy iteration

    # Define plot labels and limits
    x_label = 'Distance (Å)'
    y_label = 'g(r)'
    x_lim = (0, 6)

    for i, (atom_id, rdf_dict) in enumerate(data.items()):
        ax = axs[i]
        for label, rdf in rdf_dict.items():
            rdf_smoothed = gaussian_filter1d(rdf.rdf, sigma=2)
            ax.plot(rdf.bins, rdf_smoothed, label=label)  # Plot RDF curve

        # Title for each subplot
        ax.set_title(f"Interactions between {atom_id} and {interface_name.split('_')[0]} surface")  # Customize as needed for each RDF plot
        ax.legend()

        # Set x-axis properties, visible only for bottom row
        ax.set_xlim(x_lim)
        ax.set_ylim(bottom=0)
        if i >= (rows - 1) * cols:  # Only show x-axis labels on the last row
            ax.set_xlabel(x_label)
            ax.tick_params(axis='x', labelbottom=True)
        else:
            ax.tick_params(axis='x', labelbottom=False)

        # Set y-axis properties, visible only for leftmost plots
        if i % cols == 0:
            ax.set_ylabel(y_label)
        else:
            ax.tick_params(axis='y', labelleft=False)

    # Remove unused subplots, if any
    for j in range(num_plots, len(axs)):
        fig.delaxes(axs[j])

    # Save or display the plot
    if save_path:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(os.path.join(save_path, interface_name), bbox_inches='tight', dpi=500)
        plt.close(fig)  # Close the figure after saving
    else:
        plt.show()
def identify_element(atomsymbolfulllist, element):

    if element == None:
        element = input("What is the element that you want to identify all the ids?\n")
        while element not in atomsymbolfulllist:
            element = input("Please, provide an element that is present in the system!\n")

    ids = []
    for id, elem in enumerate(atomsymbolfulllist):
        if elem == element:
            ids.append(id)

    return ids
def analyze_hydrogen_bonds(xyz_path, x, y, z, donors_sel_ids, hydrogens_sel_ids, acceptors_sel_ids, start=1250, stop=None, step=1):
    """
    Analyzes hydrogen bonds between a specified selection of atoms within an MDAnalysis Universe.

    Args:
        xyz_path (str): Path to the XYZ file containing the molecular structure.
        x (float): Simulation box length in the x-direction.
        y (float): Simulation box length in the y-direction.
        z (float): Simulation box length in the z-direction.
        spe (str): Selection string for atoms of interest in the species (e.g., a polymer).
        surface (str): Selection string for atoms of interest in the surface (e.g., a metal oxide).
        start (int, optional): Starting frame for the analysis. Default is 1250.
        stop (int, optional): Stopping frame for the analysis. Default is None, which analyzes to the last frame.
        step (int, optional): Frame step interval for analysis. Default is 1.

    Returns:
        np.ndarray: Array of hydrogen bonds found in the specified frames. Each entry contains:
            - Donor index
            - Hydrogen index
            - Acceptor index
            - Frame index

    Notes:
        Hydrogen bonds are defined by a maximum distance of 3.5 Å and a minimum angle of 150°.
        The analysis uses MDAnalysis' `HydrogenBondAnalysis` tool to identify H-bond interactions
        between specified atom groups.
    """

    u = mda.Universe(xyz_path, format="XYZ", guess_bonds=True)
    u.dimensions = [x, y, z, 90, 90, 90]

    # Parse the XYZ file with XYZParser to infer atom types
    parser = mda.topology.XYZParser.XYZParser(xyz_path)
    topology = parser.parse()
    # print("Guessed attributes:\n-------------------")
    # for attr in topology.guessed_attributes:
    #     print(attr.attrname)
    #     if attr.attrname == 'types':
    #         print(attr.values)
    # print("Read attributes:\n-----------------")
    # for attr in topology.read_attributes:
    #     print(attr.attrname)
    #     if attr.attrname == 'names':
    #         print(attr.values)

    element_types = []
    element_masses = []
    for attr in topology.guessed_attributes:
        if attr.attrname == 'types':
            for t in attr.values:
                element_types.append(t)
        elif attr.attrname == "masses":
            for m in attr.values:
                element_masses.append(m)

    element_types = np.array(element_types)
    element_masses = np.array(element_masses)

    # Add atom types & masses to Universe
    u.add_TopologyAttr('types', element_types)
    u.add_TopologyAttr('masses', element_masses)

    print("Number of atoms :", len(u.atoms))
    print("Atoms:")
    for atom in u.atoms:
        print(atom.id, atom)
    print("Bonds found")
    for bond in u.bonds:
        print(bond)

    #donors_sel =     #"index " + " ".join(map(str, donors_sel_ids)) #AtomGroup(donors_sel_ids, u)
    #hydrogens_sel =     #"index " + " ".join(map(str, hydrogens_sel_ids)) #AtomGroup(hydrogens_sel_ids, u)
    #acceptors_sel =     #"index " + " ".join(map(str, acceptors_sel_ids)) #AtomGroup(acceptors_sel_ids, u)

    # Add bond information to the Universe
    #guess_bonds(atom, atom.positions, box=[x, y, z])
    #bonds = guess_bonds(u.atoms.positions, u.atoms.types)
    #u.add_TopologyAttr('bonds', bonds)

    # Optionally add angle information if needed for more detailed topology analysis
    #u.add_TopologyAttr('angles', guess_angles(u.bonds))

    hbonds = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donors_sel_ids,
        hydrogens_sel=hydrogens_sel_ids,
        acceptors_sel=acceptors_sel_ids,
        d_a_cutoff=3.0,
        d_h_a_angle_cutoff=150,
    )
    #hbonds.hydrogens_sel = hbonds.guess_hydrogens()
    #hbonds.acceptors_sel = hbonds.guess_acceptors()
    hbonds.run(start=0, stop=stop, step=step)
    print(f"There are {hbonds.results.hbonds.shape[0]} h-bonds in total (cumulated through frames)")
    return hbonds
def fill_gaps(hbonds):
    threshold = 25
    frame_ids = []
    total_nb_frames = len(hbonds.times)
    for i, hbond in enumerate(hbonds.results.hbonds):
        frame_ids.append(int(hbond[0]))  # Ensure frame IDs are integers

    nb_frame_with_hbonds = len(set(frame_ids))
    result = []
    for i in range(len(frame_ids) - 1):
        result.append(frame_ids[i])
        # Check the difference with the next number
        if frame_ids[i + 1] - frame_ids[i] < threshold:
            # Add the missing numbers in between, converting to int
            result.extend(range(int(frame_ids[i]) + 1, int(frame_ids[i + 1])))
    if frame_ids != []:
        result.append(frame_ids[-1])  # Add the last number
    return nb_frame_with_hbonds, result
def plot_hbond_results(hbond, save_path=None, interface_name=None):

    fig, ax = plt.subplots(1, 2, figsize=(13, 6))

    total_simulation_time = 15000  # in fs
    time_per_frame = 4 # in fs/frame

    # Custom x-tick labels
    frame_ticks = range(0, total_simulation_time, 525)
    time_ticks = [int(round(frame * time_per_frame/1000, 0)) for frame in frame_ticks]
    # Custom y-tick labels
    max_hbonds = 3
    bond_ticks = range(0, max_hbonds+2, 1)
    y_ticks = [int(round(bond/2, 0)) for bond in bond_ticks]

    ax[0].plot(hbond.times, hbond.count_by_time(), color='r')
    nb_frame_with_hbonds, filled_hbonds = fill_gaps(hbond)
    binary = np.zeros(len(hbond.times))
    for i in range(len(binary)):
        for j in filled_hbonds:
            if i==j:
                binary[i] = 1
    ax[1].plot(hbond.times, binary, color='r')

    for i in range(2):
        ax[i].set_xlabel(xlabel='Time (ps)', fontsize=18)
        ax[i].set_xticks(frame_ticks)
        ax[i].set_xticklabels(time_ticks)
        ax[i].tick_params(axis='x', labelsize=18)
        ax[i].set_xlim(left=0)
        ax[i].set_xlim(right=3750)

        ax[i].set_ylabel(ylabel="Presence of Hydrogen Bond",fontsize=18)
        ax[i].tick_params(axis='y', labelsize=18)
        ax[i].set_yticks(y_ticks)
        ax[i].set_ylim(bottom=0)
        if i == 0:
            ax[i].set_ylim(top=3)
        else:
            ax[i].set_ylim(top=1.5)


    ax[0].set_title(f"Time evolution of H-bonds number", fontsize=20)
    ax[1].set_title(f"Time evolution of H-bonds number", fontsize=20)


    fig.tight_layout()

    # Save the plot if a path and name are provided
    if save_path and interface_name:
        plt.savefig(f"{save_path}/Hbonds_{interface_name}.png", bbox_inches='tight', dpi=500)

    plt.show()
class SelectFromCollection(object):
    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other
        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))
        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []
        self.canvas.get_tk_widget().focus_set()
    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.get_tk_widget().focus_set()
        self.canvas.draw_idle()
    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

def main() -> None:
    """
    Main function to execute the XDATCAR file selection and dictionary organization.
    """

    # 1) Read all xdatcar files and put the paths into a dic
    folder_path = select_folder()
    if not folder_path:
        print("Invalid folder path. Please try again.")
        return
    xdatcar_files = get_xdatcar_files(folder_path)
    xdatcar_dict = add_paths_to_dict(xdatcar_files)

    print("Selected Folder:", folder_path)
    print("Final XDATCAR Dictionary Structure:")
    print(xdatcar_dict)

    # 2) Turn all the xdatcar files into xyz format files
    for li_status, grade in xdatcar_dict.items():
        for grade, spe in grade.items():
            for spe, xdatcar_path in spe.items():
                # Insert custom actions here
                print(f"Processing: {li_status} > {grade} > {spe}")
                print(f"File path: {xdatcar_path}")

                elect = 'NMC'
                if li_status == 'With Li':
                    elect = 'LiNMC'

                output_directory = folder_path
                xyz_file = "traj_" + elect + grade + "_" + spe
                file_name = xyz_file + ".xyz"
                output_file = os.path.join(output_directory, file_name)
                #vectorFileName = os.path.join(output_directory, xyz_file)

                lines = read_xdatcar(xdatcar_path)
                a, b, c, atomsymbolfulllist, quantities, lattice = system(lines)
                uniq_elem = list(set(atomsymbolfulllist))
                interface = write_xyz(lines, output_file, a, b, c, atomsymbolfulllist, quantities)
                #write_lattice(vectorFileName, lattice)

                # 3) Select the indices from spe and then from the electrode outer layer
                electrode_selection = manage_selection_plot('Electrode', interface, uniq_elem, radiuslist, element_colors)
                spe_selection = manage_selection_plot(spe, interface, uniq_elem, radiuslist, element_colors)

                print('Selection from electrode:', electrode_selection)
                print('Selection from SPE:', spe_selection)
                print('-'*100)

                electrode_selection_metals = []
                electrode_selection_oxygen = []
                for i in electrode_selection:
                    if atomsymbolfulllist[i] == 'O':
                        electrode_selection_oxygen.append(i)
                    else:
                        electrode_selection_metals.append(i)

                spe_selection_O_N = []
                spe_selection_hydrogen = []
                for i in spe_selection:
                    if atomsymbolfulllist[i] == 'O' or atomsymbolfulllist[i] == 'N':
                        spe_selection_O_N.append(i)
                    elif atomsymbolfulllist[i] == 'H':
                        spe_selection_hydrogen.append(i)

                data = {}
                # 4) Iterate through all spe_ind / electrode_ind comb and calculate rdf
                for ind_spe in spe_selection_O_N:
                    for ind_electr in electrode_selection_metals:
                        spe_id = atomsymbolfulllist[ind_spe] + str(ind_spe)   # e.g. C2
                        metal_id = atomsymbolfulllist[ind_electr] + str(ind_electr)   # e.g. Li25
                        label = spe_id + "-" + metal_id

                        rdf = rdf_calculation(output_file, a, b, c, ind_spe, ind_electr)

                        if not rdf_meets_criterion(rdf, criterion_range=3):
                            if spe_id not in data and li_status == 'With Li':
                                data[spe_id] = {}  # Initialize as an empty dictionary if spe_id is not in data
                            #rdf_smoothed = gaussian_filter1d(rdf, sigma=2)
                            data[spe_id][label] = rdf#_smoothed
                            print("Data to plot:")
                            print(data)

                # 5) Plot rdf for the function that are non-0 before 3
                plot_rdfs(data, save_path=output_directory, interface_name=elect+grade+"_"+spe+"_smoothened")

                # 6) Calculate H-bonds count and plot its time evolution (only for N/O-H from Alcohol, Amine and Urethane)
                # if spe == 'Alcohol':
                #     donors_sel_ids = f"index {len(atomsymbolfulllist)-3}:{len(atomsymbolfulllist)-1}"
                #     hydrogens_sel_ids = f"index {identify_element(atomsymbolfulllist, 'H')[0]}:{identify_element(atomsymbolfulllist, 'H')[-1]}"
                #     acceptors_sel_ids = f"index {identify_element(atomsymbolfulllist, 'O')[0]}:{identify_element(atomsymbolfulllist, 'O')[-4]}"
                # elif spe == 'Amine':
                #     donors_sel_ids = f"index {identify_element(atomsymbolfulllist, 'N')[0]}:{identify_element(atomsymbolfulllist, 'N')[-1]}"
                #     hydrogens_sel_ids = f"index {identify_element(atomsymbolfulllist, 'H')[-2]}:{identify_element(atomsymbolfulllist, 'H')[-1]}"
                #     acceptors_sel_ids = f"index {identify_element(atomsymbolfulllist, 'O')[0]}:{identify_element(atomsymbolfulllist, 'O')[-1]}"
                # elif spe == 'Urethane':
                #     donors_sel_ids = f"index {identify_element(atomsymbolfulllist, 'N')[0]}:{identify_element(atomsymbolfulllist, 'N')[-1]}"
                #     hydrogens_sel_ids = f"index {identify_element(atomsymbolfulllist, 'H')[-2]}:{identify_element(atomsymbolfulllist, 'H')[-1]}"
                #     acceptors_sel_ids = f"index {identify_element(atomsymbolfulllist, 'O')[0]}:{identify_element(atomsymbolfulllist, 'O')[-5]}"
                #
                # if spe in ['Alcohol', 'Amine', 'Urethane']:
                #     hbond_data = analyze_hydrogen_bonds(output_file, a, b, c, donors_sel_ids, hydrogens_sel_ids, acceptors_sel_ids)
                #     plot_hbond_results(hbond_data, save_path=output_directory, interface_name=elect+grade+"_"+spe)


if __name__ == "__main__":
    main()
