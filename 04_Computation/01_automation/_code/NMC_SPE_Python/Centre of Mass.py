import mplcursors
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Button
from tkinter import filedialog
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.text import Annotation

# Global variable for the axes
ax = None

def read_input_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def parse_atoms_info(lines):
    scaling_factor = float(lines[1])
    cell_dimensions = np.array([list(map(float, line.split())) for line in lines[2:5]])
    elements = lines[5].split()
    print("Cell dimmensions : ", cell_dimensions)
    return elements, scaling_factor, cell_dimensions

def parse_atoms_count(lines):
    atoms = lines[5].split()
    counts = list(map(int, lines[6].split()))
    print(counts)
    return dict(zip(atoms, counts))

def parse_atom_coordinates(lines):
    start_line = 9  # Skipping unnecessary lines
    coordinates = np.array([list(map(float, line.split()[:3])) for line in lines[start_line:]])
    return coordinates

def calculate_center_of_mass(elements, scaling_factor, cell_dimensions, counts, coordinates):
    molarweightlist = {"H": 1.00794, 'Li': 6.941, "O": 15.9994, "C": 12.0107, "Ar": 39.948, "He": 4.002602, "Ne": 20.1797,
                       "N": 14.0067, "Cl": 35.453, "Na": 22.989769, "F": 18.998403, "Br": 79.904, "Tc": 98.0,
                       "Pm": 145.0, "Cd": 112.411, "Mn": 54.9380, "Ni": 58.70, "Co": 58.9332, "fulvic_acid": 308.242, "Pb": 207.2, 'Cu': 63.546,
                       "Fe": 55.845}

    total_mass = 0
    center_of_mass = np.zeros(3)

    cumulative_counts = np.cumsum(list(counts.values()))
    print("Cumulative counts : ", cumulative_counts)
    a, b, c = cell_dimensions

    for element, count in counts.items():
        mass = molarweightlist[element]
        total_mass += mass * count

        if element in elements:
            print(element)
            element_indices = np.arange(cumulative_counts[elements.index(element)] - count, cumulative_counts[elements.index(element)])
            #print(element_indices)
            element_coordinates = coordinates[element_indices]
            #print(element_coordinates)
            center_of_mass += np.sum(element_coordinates * mass, axis=0)
            print(center_of_mass)

    center_of_mass /= total_mass
    center_of_mass = center_of_mass * scaling_factor  # Scaling back to original dimensions

    return center_of_mass

def extract_atom_data(lines):
    start_line = 9  # Skipping unnecessary lines
    atom_data = []

    # Parse the element names and occurrences
    element_names = lines[5].split()
    element_occurrences = list(map(int, lines[6].split()))

    # Combine element names and their occurrences into a single list
    elements = [element for name, count in zip(element_names, element_occurrences) for element in [name] * count]

    for line, element in zip(lines[start_line:], elements):
        elements_and_coordinates = line.split()
        # Extracting all coordinates (including negative ones)
        coordinates = [float(coord) for coord in elements_and_coordinates[:3]]
        atom_data.append([element] + coordinates)

    return atom_data


def rotate_x_pos(event):
    ax.view_init(elev=ax.elev, azim=ax.azim - 45)
    plt.draw()
def rotate_x_neg(event):
    ax.view_init(elev=ax.elev, azim=ax.azim + 45)
    plt.draw()
def rotate_y_pos(event):
    ax.view_init(elev=ax.elev - 45, azim=ax.azim)
    plt.draw()
def rotate_y_neg(event):
    ax.view_init(elev=ax.elev + 45, azim=ax.azim)
    plt.draw()
def rotate_z_pos(event):
    ax.view_init(elev=ax.elev, azim=ax.azim - 45)
    plt.draw()
def rotate_z_neg(event):
    ax.view_init(elev=ax.elev, azim=ax.azim + 45)
    plt.draw()
def align_view_x(event):
    global ax
    ax.view_init(elev=0, azim=0)
    plt.draw()
def align_view_y(event):
    global ax
    ax.view_init(elev=0, azim=-90)
    plt.draw()
def align_view_z(event):
    global ax
    ax.view_init(elev=90, azim=0)
    plt.draw()
def visualize_atoms_and_com(atom_data, element_colors, radiuslist):
    global ax
    fig = plt.figure(figsize=(10, 5))

    # Main 3D scatter plot
    ax = fig.add_subplot(121, projection='3d')

    # Lists to store legend handles and labels
    marker_handles = []
    marker_labels = []

    # Plotting atoms
    for entry in atom_data:
        element = entry[0]
        coordinates = entry[1:]
        color = element_colors.get(element, 'black')  # Use 'black' if color not defined
        radius = radiuslist.get(element, 1.0)  # Use 1.0 if radius not defined
        ax.scatter(coordinates[0], coordinates[1], coordinates[2], label=element, c=color, marker='o', s=radius**2 * 100)

        # Add legend handles and labels for markers
        if element not in marker_labels:
            marker_handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                                         markersize=8, label=element))  # Adjust the markersize value as needed
            marker_labels.append(element)

        # Create custom legends for materials and markers
        marker_legend = plt.legend(handles=marker_handles, labels=marker_labels, title='Element', loc='upper left',
                                       bbox_to_anchor=(0, 1), fontsize=10, markerscale=1.5, title_fontsize=12)

    # Plotting center of mass
    ax.scatter(center_of_mass[0], center_of_mass[1], center_of_mass[2], c='r', marker='x', s=100, label='Center of Mass')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Creating legend with only the elements present in the system
    ax.add_artist(marker_legend)

    # Add buttons for aligning the view
    mplcursors.cursor(hover=True)
    ax.annotate('Align X', xy=(1, 0.5), xytext=(1.1, 0.5), arrowprops=dict(facecolor='black', shrink=0.05))
    ax.annotate('Align Y', xy=(0.5, 1), xytext=(0.5, 1.1), arrowprops=dict(facecolor='black', shrink=0.05))
    ax.annotate('Align Z', xy=(0.5, 0), xytext=(0.5, -0.1), arrowprops=dict(facecolor='black', shrink=0.05))

    # Add buttons for aligning the view
    ax_align_x = plt.axes([0.65, 0.01, 0.1, 0.05])
    ax_align_y = plt.axes([0.76, 0.01, 0.1, 0.05])
    ax_align_z = plt.axes([0.87, 0.01, 0.1, 0.05])
    btn_align_x = Button(ax_align_x, 'Align X')
    btn_align_y = Button(ax_align_y, 'Align Y')
    btn_align_z = Button(ax_align_z, 'Align Z')
    btn_align_x.on_clicked(align_view_x)
    btn_align_y.on_clicked(align_view_y)
    btn_align_z.on_clicked(align_view_z)

    # Add buttons for rotating around each axis
    ax_rotate_x_pos = plt.axes([0.05, 0.06, 0.1, 0.04])
    ax_rotate_x_neg = plt.axes([0.05, 0.01, 0.1, 0.04])
    ax_rotate_y_pos = plt.axes([0.2, 0.06, 0.1, 0.04])
    ax_rotate_y_neg = plt.axes([0.2, 0.01, 0.1, 0.04])
    ax_rotate_z_pos = plt.axes([0.35, 0.06, 0.1, 0.04])
    ax_rotate_z_neg = plt.axes([0.35, 0.01, 0.1, 0.04])
    btn_rotate_x_pos = Button(ax_rotate_x_pos, '+45°')
    btn_rotate_x_neg = Button(ax_rotate_x_neg, '-45°')
    btn_rotate_y_pos = Button(ax_rotate_y_pos, '+45°')
    btn_rotate_y_neg = Button(ax_rotate_y_neg, '-45°')
    btn_rotate_z_pos = Button(ax_rotate_z_pos, '+45°')
    btn_rotate_z_neg = Button(ax_rotate_z_neg, '-45°')
    btn_rotate_x_pos.on_clicked(rotate_x_pos)
    btn_rotate_x_neg.on_clicked(rotate_x_neg)
    btn_rotate_y_pos.on_clicked(rotate_y_pos)
    btn_rotate_y_neg.on_clicked(rotate_y_neg)
    btn_rotate_z_pos.on_clicked(rotate_z_pos)
    btn_rotate_z_neg.on_clicked(rotate_z_neg)

    # Add text labels for each column
    ax_text_labels = fig.add_axes([0.0, 0.1, 0.8, 0.05], frame_on=False)
    ax_text_labels.text(0.13, 0.3, 'X-Axis', fontsize=10, ha='center')
    ax_text_labels.text(0.32, 0.3, 'Y-Axis', fontsize=10, ha='center')
    ax_text_labels.text(0.50, 0.3, 'Z-Axis', fontsize=10, ha='center')
    ax_text_labels.set_xticks([])
    ax_text_labels.set_yticks([])


    # Small subplot for graphical representation of axes
    ax_representation = fig.add_subplot(122, projection='3d')
    ax_representation.set_box_aspect([1, 1, 1])  # Set the aspect ratio to be equal
    ax_representation_position = [0.5, 0.1, 0.2, 0.8]  # [left, bottom, width, height]
    ax_representation.set_position(ax_representation_position)
    ax_representation.set_axis_off()

    # Draw visual representation of axes in the small subplot
    arrow_length = 0.03
    arrow_size = 0.01

    # X-axis
    ax_representation.quiver(0, 0, 0, arrow_length, 0, 0, color='black', arrow_length_ratio=0.1)
    ax_representation.text(arrow_length + arrow_size, 0, 0, 'X', color='black', ha='center', va='center', fontsize=10)

    # Y-axis
    ax_representation.quiver(0, 0, 0, 0, arrow_length, 0, color='black', arrow_length_ratio=0.1)
    ax_representation.text(0, arrow_length + arrow_size, 0, 'Y', color='black', ha='center', va='center', fontsize=10)

    # Z-axis
    ax_representation.quiver(0, 0, 0, 0, 0, arrow_length, color='black', arrow_length_ratio=0.1)
    ax_representation.text(0, 0, arrow_length + arrow_size, 'Z', color='black', ha='center', va='center', fontsize=10)

    # Synchronize the view in both subplots
    def update_representation_view(event):
        ax_representation.view_init(elev=ax.elev, azim=ax.azim)
        fig.canvas.draw()

    fig.canvas.mpl_connect('motion_notify_event', update_representation_view)
    plt.show()



if __name__ == "__main__":
    # Define a dictionary of colors for each element
    element_colors = {'H': 'white', 'Li': 'green', 'O': 'red', 'C': 'brown', 'N': 'blue',
                      'Mn': 'pink', 'Ni': 'grey', 'Co': 'indigo'}

    radiuslist = {"H": 1.20, "Li": 1.82, "Cl": 2.75, "C": 1.70, "Si": 2.10, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80,
                  'Ne': 1.54, 'Ar': 1.88, 'He': 1.40, "Na": 2.27, "F": 1.47, "Br": 1.85, "Tc": 2.05, "Pm": 2.43,
                  "Cd": 1.55, "Co": 1.35, "Mn": 1.61, "Fe": 1.61,
                  "Ni": 1.49}  # IN ANGSTROM (Molecular coord also in A) (from http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html)

    file_path = filedialog.askopenfilename(title="Select the POSCAR/CONTCAR file")
    lines = read_input_file(file_path)

    elements, scaling_factor, cell_dimensions = parse_atoms_info(lines)
    counts = parse_atoms_count(lines)
    print("Elements:", elements)
    print("counts : ", counts)
    coordinates = parse_atom_coordinates(lines)

    center_of_mass = calculate_center_of_mass(elements, scaling_factor, cell_dimensions, counts, coordinates)

    print("Center of Mass:", center_of_mass)

    atom_data = extract_atom_data(lines)
    visualize_atoms_and_com(atom_data, element_colors, radiuslist)