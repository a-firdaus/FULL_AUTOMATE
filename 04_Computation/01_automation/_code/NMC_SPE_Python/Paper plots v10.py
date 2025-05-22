import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
import matplotlib.colors as mcolors

# This v10 plots data 2 datasets in the same graph but with different shade
# (to compare with and w/o Li or with and w/o EF)

def custom_sort(data):
    # Define the sorting key function
    def sorting_key(item):
        # Custom sorting order
        sorting_order = {'C-C':0, 'C-H':1, 'C-N (amine)':2, 'C-N (urethane)':3,
                         'C-O (carbonate)':4, 'C-O (ester)':5, 'C-O (ether)':6, 'C-O (urethane)':7,
                         'C-O (alcohol)':8, 'C=O (carbonate)':9, 'C=O (ester)':10, 'C=O (urethane)':11,
                         'N-H (amine)':12, 'N-H (urethane)':13, 'O-H':14}
        # Use the sorting order value as the key for sorting
        return sorting_order.get(item[0], float('inf'))

    # Sort the data using the custom sorting key
    sorted_data = sorted(data, key=sorting_key)

    return sorted_data
def sort_Mode(data):
    # Define the sorting key function
    def sorting_key(item):
        # Custom sorting order
        sorting_order = {'C-C':0, 'C-H':1, 'C-N (amine)':10, 'C-N (urethane)':8,
                         'C-O (carbonate)':4, 'C-O (ester)':12, 'C-O (ether)':14, 'C-O (urethane)':6,
                         'C-O (alcohol)':3, 'C=O (carbonate)':5, 'C=O (ester)':13, 'C=O (urethane)':7,
                         'N-H (amine)':11, 'N-H (urethane)':9, 'O-H':2}
        # Use the sorting order value as the key for sorting
        return sorting_order.get(item[0], float('inf'))

    # Sort the data using the custom sorting key
    sorted_data = sorted(data, key=sorting_key)

    return sorted_data
def input_name_file():
    nameFile = input('Name of the file to save: ')
    return nameFile
def input_comp():
    comp = ''
    print('Do you want to plot some data onto the same plot with shade of colours?')
    while comp.upper() != 'YES' and comp.upper() != 'NO':
        comp = input()
    return comp
def input_several_data():
    print('How many data do you want to plot?')
    while True:
        try:
            several_data = int(input())
            # If the user enters a non-integer value, the ValueError will be raised
            break  # Break out of the loop if input is successful
        except ValueError:
            print("Invalid input. Please enter a valid integer.")
    return several_data
def input_data():
    data = ''
    print("Do you want to plot data with (Enter 'with') or without (Enter 'without') electric field?")
    while data.upper() != 'WITH' and data.upper() != 'WITHOUT':
        data = input()
    return data
def input_comp_type():
    print(
        "What comparison do you want to make? Presence of Li (Enter 'Li'); Influence of Electric field (Enter 'EF').")
    comp_type = input()
    while comp_type.upper() not in ['LI', 'EF']:
        print('Invalid input. Please enter a new input.')
        comp_type = input()
    return comp_type
def append_multiple_entries(entries, final_list):
    for entry in entries:
        final_list.append(entry)
    return final_list
def append_materials(comp, comp_type, efOrNot, several_data, names_data):

    materials = []

    if comp.upper() == 'YES':
        if comp_type.upper() == 'LI':
            #append_multiple_entries([fourth_data, third_data, data, data_bis], materials)
            append_multiple_entries(['NMC622', 'LiNMC622', 'NMC811', 'LiNMC811'], materials)

        else:
            #ppend_multiple_entries(['NMC622', 'NMC622 EF', 'NMC811', 'NMC811 EF'], materials)
            append_multiple_entries(['LiNMC622', 'LiNMC622 EF', 'LiNMC811', 'LiNMC811 EF'], materials)

    else:
        if efOrNot.upper() == 'WITH':
            for i in range(several_data):
                while True:
                    try:
                        new_mat = input(f"Enter material {i + 1}: ") + ' EF'
                        print(new_mat)
                        if new_mat not in names_data:
                            raise ValueError("Material not defined in the dictionary.")
                        if new_mat in materials:
                            raise ValueError("Material already selected.")

                        materials.append(new_mat)
                        break  # Break out of the loop if input is successful
                    except ValueError as e:
                        print(f"Error: {e} Please enter a valid material.")
        else:
            for j in range(several_data):
                while True:
                    try:
                        new_mat = input(f"Enter material {j + 1}: ")
                        if new_mat not in names_data:
                            raise ValueError("Material not defined in the dictionary.")
                        if new_mat in materials:
                            raise ValueError("Material already selected.")

                        materials.append(new_mat)
                        break  # Break out of the loop if input is successful
                    except ValueError as e:
                        print(f"Error: {e} Please enter a valid material.")

    print("Selected materials:", materials)
    return materials
def append_mat_dico(materials, names_data):
    material_data_dico = {}
    for mat in materials:
        names_data_for_material = names_data[mat]

        # Create variables dynamically
        names = [item[0] for item in names_data_for_material]
        min_values = [item[3] for item in names_data_for_material]
        max_values = [item[2] for item in names_data_for_material]
        mode_values = [item[1] for item in names_data_for_material]

        # Store the lists in the material_data dictionary :
        # dico whose keys are mat (e.g. 'LiNMC622 EF') and values are dico whose keys are what can be found
        # in data (i.e. bond name, min, max and mode) and values are list of values
        material_data_dico[mat] = {
            'names': names,
            'min_values': min_values,
            'max_values': max_values,
            'mode_values': mode_values
        }

        print(material_data_dico[mat])


    return material_data_dico
def sort_data(several_data, materials, names_data):

    sorted_names_list = []
    sorted_min_values_list = []
    sorted_max_values_list = []
    sorted_mode_values_list = []

    for i in range(several_data):
        # Extracting data for the three plots using sorted_data
        sorted_data = custom_sort(names_data[materials[i]])

        # Append the sorted names and values to the respective lists
        sorted_names_list.append([item[0] for item in sorted_data])
        sorted_min_values_list.append([item[3] for item in sorted_data])
        sorted_max_values_list.append([item[2] for item in sorted_data])
        sorted_mode_values_list.append([item[1] for item in sorted_data])


    # Create lists for sorted names and values dynamically
    names_list = []
    min_values_list = []
    max_values_list = []
    mode_values_list = []

    # Create a set to keep track of unique names while preserving order
    ordered_set = []
    seen_names = set()


    for i in range(several_data):

        # Extracting data for the three plots using sorted_data
        data = names_data[materials[i]]

        # Append the sorted names and values to the respective lists
        names_list.append([item[0] for item in data])
        min_values_list.append([item[3] for item in data])
        max_values_list.append([item[2] for item in data])
        mode_values_list.append([item[1] for item in data])

    # Iterate over sorted_data to maintain order
    for item in sorted_data:
        name = item[0]
        if name not in seen_names:
            ordered_set.append(name)
            seen_names.add(name)

    return ordered_set, sorted_names_list, sorted_min_values_list, sorted_max_values_list, sorted_mode_values_list
def draw_several_subplots(materials, names_data):

    colors ={'LiNMC622': mcolors.TABLEAU_COLORS['tab:orange'],
             'LiNMC811': mcolors.TABLEAU_COLORS['tab:blue'],
             'NMC811': mcolors.TABLEAU_COLORS['tab:red'],
             'NMC622': mcolors.TABLEAU_COLORS['tab:green'],
             'LiNMC622 EF': mcolors.TABLEAU_COLORS['tab:orange'],
             'LiNMC811 EF': mcolors.TABLEAU_COLORS['tab:blue'],
             'NMC811 EF': mcolors.TABLEAU_COLORS['tab:red'],
             'NMC622 EF': mcolors.TABLEAU_COLORS['tab:green']
             }

    num_columns = min(2, len(materials))  # Set num_columns to 2 or the length of materials, whichever is smaller
    num_rows = (len(materials) + num_columns - 1) // num_columns

    fig, axs = plt.subplots(num_rows, num_columns, figsize=(12, 6), sharex=True, sharey=True, constrained_layout=True)

    # Ensure axs is a 2D array
    if num_rows == 1:
        axs = axs.reshape(1, -1)
    elif num_columns == 1:
        axs = axs.reshape(-1, 1)

    for i, material in enumerate(materials):
        if i >= num_columns * num_rows:  # Only iterate up to the number of subplots
            break

        row_index = i // num_columns
        col_index = i % num_columns
        print(row_index, col_index)

        names_data_for_material = names_data[material]
        x_data = [item[3] for item in names_data_for_material]
        y_data = [item[2] for item in names_data_for_material]

        axs[row_index, col_index].scatter(x_data, y_data, color=colors[material.split()[0]], marker='o')
        axs[row_index, col_index].set_xlabel('MinRBC (%)', fontsize=12, fontweight='bold')
        if col_index == 0:
            axs[row_index, col_index].set_ylabel('MaxRBC (%)', fontsize=12, fontweight='bold')
        axs[row_index, col_index].set_title(material)
        axs[row_index, col_index].set_ylim(bottom=0)
        axs[row_index, col_index].set_xlim(left=-14)
        axs[row_index, col_index].set_xlim(right=0)
        axs[row_index, col_index].set_ylim(top=30)
        axs[row_index, col_index].yaxis.labelpad = 10
        axs[row_index, col_index].xaxis.labelpad = 10
        axs[row_index, col_index].grid(True, linestyle='--', alpha=0.7)


    plt.tight_layout()
    print('Save the file (Max vs Min)\n--------------------------')
    plt.savefig(input_name_file() + '.png', bbox_inches='tight', dpi=500)
    #plt.show()

    # Create an empty DataFrame to store the data
    df_whisker = pd.DataFrame()

    # Add data for each material to the DataFrame
    for i, mat in enumerate(materials):
        df_material = pd.DataFrame({
            'Bond type': material_data_dico[mat]['names'],
            'Mode Values': material_data_dico[mat]['mode_values'],
            'Material': [mat] * len(material_data_dico[mat]['names'])
        })
        df_whisker = pd.concat([df_whisker, df_material], ignore_index=True)

    # Plotting the whisker-type chart using seaborn
    fig, ax = plt.subplots(figsize=(10, 5))
    # Use the material_colors dictionary for setting custom colors
    sns.boxplot(x='Bond type', y='Mode Values', hue='Material', data=df_whisker, ax=ax,
                palette=colors, whis=[0, 100])
    ax.set_xlabel('Bond type', fontweight='bold', fontsize=16)
    ax.set_ylabel('ModeRBC (%)', fontweight='bold', fontsize=16)

    ax.set_ylim(bottom=-5)
    ax.set_ylim(top= 5)

    ax.set_xticks(df_whisker['Bond type'])
    ax.set_xticklabels(df_whisker['Bond type'])
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    # Adjust the space between axis and their label
    ax.yaxis.labelpad = 14  # Adjust this value as needed
    ax.xaxis.labelpad = 14  # Adjust this value as needed

    # Adjust layout
    plt.tight_layout()
    print('Save the file (Mode)\n--------------------------')
    plt.savefig(input_name_file() + '.png', bbox_inches='tight', dpi=500)
    plt.show()
def draw_subplot_single_mat(mat, material_data_dico, ordered_set, sorted_names_list, sorted_min_values_list, sorted_max_values_list):
    y_max = 30
    # Set the size of the points
    point_size = 10
    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(4, 6, width_ratios=[1, 1, 1, 1, 0.2, 0.2], height_ratios=[1, 1, 1, 1])

    # Iterate over unique names
    for i, unique_name in enumerate(ordered_set):
        exceed_limit = False
        # Filter data for the current unique name using sorted_data
        filtered_min = [min_val for names, min_vals in zip(sorted_names_list, sorted_min_values_list) for name, min_val
                        in zip(names, min_vals) if name == unique_name]
        filtered_max = [max_val for names, max_vals in zip(sorted_names_list, sorted_max_values_list) for name, max_val
                        in zip(names, max_vals) if name == unique_name]

        # Check if any data point exceeds y_max
        if any(value > y_max for value in filtered_max):
            exceed_limit = True

        # Calculate grid position
        row = i // 4
        col = i % 4

        # Plot the data on the current subplot with the same color
        ax = plt.subplot(gs[row, col])
        ax.scatter(filtered_min, filtered_max, s=point_size, c='orange', marker='o')

        # Share y-axis label for each row
        if col == 0:
            ax.set_ylabel('MaxRBC (%)', fontsize=8, fontweight='bold', labelpad=10)

        # Share x-axis label for each column
        if row == 3 or (row == 2 and col == 3):
            ax.set_xlabel('MinRBC (%)', fontsize=8, fontweight='bold')

        ax.set_title(unique_name, fontsize=8)

        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)

        # Set y-axis limits to include zero
        ax.set_ylim(bottom=0)
        # Set y-axis limits to limit the MaxRBC displayed on the plot
        ax.set_ylim(top=y_max)
        # Set x-axis limits to include zero
        ax.set_xlim(right=0)
        # Set x-axis limits to the desired value to get the same x-axis scale for every subplot
        ax.set_xlim(left=-14)

        # If any data point exceeds y_max, mark subplot with a red asterisk
        if exceed_limit:
            ax.text(0.95, 0.90, '*', color='red', fontsize=16, fontweight='bold', ha='center', va='center',
                    transform=ax.transAxes)

    # Hide the empty subplots
    for i in range(4, 6):
        plt.subplot(gs[:, i]).set_visible(False)

    # Adjust layout
    plt.tight_layout()
    print('Save the file (Max vs Min)\n--------------------------')
    plt.savefig(input_name_file() + '.png', bbox_inches='tight', dpi=500)
    # Show the plot
    plt.show()

    names = material_data_dico[mat]['names']
    mode_values = material_data_dico[mat]['mode_values']

    # Plotting the whisker-type chart using seaborn
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(x=names, y=mode_values, ax=ax, whis=[0, 100])
    ax.set_xlabel('Bond type', fontweight='bold')
    ax.set_ylabel('ModeRBC (%)', fontweight='bold')
    #ax.set_title(uniqueNames)

    ax.set_xticks(names)
    ax.set_xticklabels(names)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    # Adjust layout
    plt.tight_layout()
    print('Save the file (Mode)\n--------------------------')
    plt.savefig(input_name_file() + '.png', bbox_inches='tight', dpi=500)
    plt.show()
def draw_comparative_subplots(comp_type, material_data_dico):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    y_max = 20

    if comp_type.upper() == 'LI':
        mat_sets = [
            (['NMC622', 'LiNMC622'], axes[0], '622', mcolors.TABLEAU_COLORS['tab:orange']),
            (['NMC811', 'LiNMC811'], axes[1], '811', 'blue')  #mcolors.TABLEAU_COLORS['tab:blue']
        ]
    elif comp_type.upper() == 'EF':
         mat_sets = [
             (['LiNMC622', 'LiNMC622 EF'], axes[0], '622', mcolors.TABLEAU_COLORS['tab:orange']),
             (['LiNMC811', 'LiNMC811 EF'], axes[1], '811', mcolors.TABLEAU_COLORS['tab:blue'])
         ]
        #mat_sets = [
        #         (['NMC622', 'NMC622 EF'], axes[0], '622', mcolors.TABLEAU_COLORS['tab:orange']),
        #         (['NMC811', 'NMC811 EF'], axes[1], '811', mcolors.TABLEAU_COLORS['tab:blue'])
        #     ]
    else:
        raise ValueError("Unknown comparison type. Please use 'LI' or 'EF'.")

    for mats, ax, title, color_base in mat_sets:
        exceed_limit = False  # Flag to check if any data point exceeds y_max
        for i, mat in enumerate(mats):
            if mat in material_data_dico:
                data = material_data_dico[mat]
                numb_data_set = len(data['names'])
                print(numb_data_set, 'Data points in', mat)
                shade = [0.3 + i * 0.5]*numb_data_set  # Different shades for different materials
                ax.scatter(data['min_values'], data['max_values'], alpha=shade, label=mat, color=color_base)

                # Check if any data point exceeds y_max
                if any(value > y_max for value in data['max_values']):
                    exceed_limit = True

            axes[i].tick_params(axis='both', which='major', labelsize=16)


        ax.set_ylim(bottom= 0,top=y_max)
        ax.set_xlim(left=-14, right=0)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.set_title(f"Material {title}")
        ax.set_xlabel('MinRBC (%)', fontsize=20)
        ax.legend(loc='upper left', fontsize=18)

        # If any data point exceeds y_max, mark subplot with a red asterisk
        if exceed_limit:
            ax.text(0.95, 0.95, '*', color='red', fontsize=20, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

    axes[1].set_yticklabels([])
    axes[0].set_ylabel('MaxRBC (%)', fontsize=20)

    print('Save the file (Max vs Min)\n--------------------------')
    plt.savefig(input_name_file() + '.png', bbox_inches='tight', dpi=500)
    plt.tight_layout()
    plt.show()

def plot_jointplots(material_data_dico, comp_type):
    if comp_type.upper() == 'LI':
        mat_sets = [
            (['NMC622', 'LiNMC622'], '622', 'Oranges'),
            (['NMC811', 'LiNMC811'], '811', 'Blues')
        ]
    elif comp_type.upper() == 'EF':
        mat_sets = [
            (['LiNMC622', 'LiNMC622 EF'], '622', 'Oranges'),
            (['LiNMC811', 'LiNMC811 EF'], '811', 'Blues')
        ]
    else:
        raise ValueError("Unknown comparison type. Please use 'LI' or 'EF'.")

    for mats, title, color in mat_sets:
        data = []
        for mat in mats:
            if mat in material_data_dico:
                df = pd.DataFrame({
                    'Min Values': material_data_dico[mat]['min_values'],
                    'Max Values': material_data_dico[mat]['max_values'],
                    'Material': mat
                })
                data.append(df)

        if data:
            combined_df = pd.concat(data, ignore_index=True)

            jointplot = sns.jointplot(data=combined_df, x='Min Values', y='Max Values', hue='Material', kind='scatter', palette=color)
            jointplot.fig.tight_layout(rect=[0, 0, 1, 1])

            print('Save the file (Mode)\n--------------------------')
            plt.savefig(input_name_file() + '.png', bbox_inches='tight', dpi=500)
            # Save or show the plot
            plt.show()

if __name__ == "__main__":

    data = sort_Mode([
        ['C-H', 1.44, 7.1, -2.58],
        ['C-H', 2.12, 6.35, -2.02],
        ['C-O (carbonate)', -1.05, 6.12, -5.94],
        ['C-O (carbonate)', -1.26, 7.21, -7.79],
        ['C-H', -0.84, 1.99, -3.5],
        ['C-H', 0.67, 4.34, -1.62],
        ['C-H', -0.85, 1.59, -3.78],
        ['C-O (carbonate)', -0.84, 6.41, -6.43],
        ['C-H', 0.71, 4.66, -2.73],
        ['C-H', 0.84, 5.28, -2.45],
        ['C-H', 0.47, 5.85, -3.2],
        ['C-O (carbonate)', -2.81, 4.85, -7.61],
        ['C-O (carbonate)', -1.11, 3.74, -5.49],
        ['C-O (carbonate)', -1.55, 6.11, -6.55],
        ['C=O (carbonate)', -1.15, 2.06, -3.52],
        ['C-O (carbonate)', -0.82, 6.96, -5.9],
        ['C-O (carbonate)', -1.2, 4.79, -5.28],
        ['C=O (carbonate)', -0.43, 1.92, -2.91],
        ['C-H', 0.7, 4.43, -3.11],
        ['C-H', 0.29, 4.84, -2.94],
        ['C-H', 0.23, 5.11, -3.23],
        ['C-O (urethane)', -1.6, 5.8, -6.84],
        ['C-H', 1.5, 7.06, -3.54],
        ['C-H', 0.69, 5.17, -3.35],
        ['C-N (urethane)', -0.64, 5.21, -4.95],
        ['C-O (urethane)', -1.07, 7.22, -7.99],
        ['C-N (urethane)', 0.78, 5.11, -4.01],
        ['C-O (urethane)', -1.94, 3.61, -6.76],
        ['C=O (urethane)', -0.02, 3.43, -3.4],
        ['C-N (urethane)', -0.62, 5.36, -4.37],
        ['C-O (urethane)', -1.65, 4.34, -6.65],
        ['C=O (urethane)', 0.1, 5.18, -2.81],
        ['C-H', 0.03, 4.66, -3.13],
        ['C-H', 0.3, 4.54, -2.97],
        ['C-H', -0.14, 4.75, -3.48],
        ['C-N (urethane)', -0.32, 5.39, -5.09],
        ['N-H (urethane)', 0.22, 3.77, -2.63],
        ['N-H (urethane)', 0.44, 6.65, -2.69],
        ['C-H', -0.59, 4.8, -5.09],
        ['C-H', -2.22, 3.9, -6.22],
        ['C-H', 1.12, 7.78, -3.41],
        ['C-N (amine)', -1.15, 5.39, -7.82],
        ['C-H', -2.02, 8.92, -7.53],
        ['C-H', 1.86, 17.31, -5.25],
        ['C-N (amine)', -1.37, 5.39, -7.24],
        ['C-N (amine)', -1.09, 7.59, -7.5],
        ['C-H', 2.81, 12.34, -3.68],
        ['C-H', 0.58, 9.68, -5.12],
        ['C-N (amine)', -1.88, 6.45, -6.67],
        ['C-N (amine)', -1.12, 4.99, -6.43],
        ['C-H', 0.34, 6.37, -4.3],
        ['C-H', 1.88, 7.22, -3.58],
        ['C-H', -0.16, 5.15, -4.6],
        ['C-N (amine)', -0.66, 6.63, -6.27],
        ['N-H (amine)', 0.06, 5.48, -4.46],
        ['N-H (amine)', 0.18, 4.99, -2.7],
        ['N-H (amine)', -0.62, 5.25, -4.96],
        ['C-H', -0.16, 4.3, -4.04],
        ['C-H', -0.92, 4.31, -5.12],
        ['C-H', 1.12, 5.56, -2.59],
        ['C-O (ester)', -2.43, 3.98, -7.37],
        ['C-C', 0.23, 7.67, -4.64],
        ['C-H', -1.26, 3.83, -4.44],
        ['C-H', 1.69, 6.27, -2.14],
        ['C-O (ester)', -0.63, 5.31, -6.0],
        ['C-O (ester)', -2.44, 4.9, -6.14],
        ['C=O (ester)', 0.26, 2.8, -2.66],
        ['C-C', 0.36, 7.29, -6.17],
        ['C-H', 1.03, 6.54, -3.22],
        ['C-H', -1.06, 4.1, -6.12],
        ['C-O (ester)', -1.47, 4.81, -6.92],
        ['C-O (ester)', -1.44, 6.75, -6.79],
        ['C=O (ester)', -0.16, 2.85, -2.38],
        ['C-C', 0.0, 5.77, -6.39],
        ['C-H', 0.23, 5.38, -4.44],
        ['C-H', 2.7, 8.77, -1.6],
        ['C-H', -0.01, 4.78, -4.35],
        ['C-O (ester)', -0.63, 6.38, -6.24],
        ['C=O (ester)', -0.26, 3.97, -3.27],
        ['C-H', 1.41, 7.82, -3.33],
        ['C-H', 0.57, 6.96, -3.79],
        ['C-H', 0.44, 5.49, -3.44],
        ['C-O (ether)', -0.76, 5.95, -6.84],
        ['C-H', 0.05, 6.93, -5.24],
        ['C-H', -0.04, 7.02, -5.48],
        ['C-O (ether)', 0.12, 7.33, -6.4],
        ['C-O (ether)', -0.4, 7.5, -7.84],
        ['C-H', -1.02, 5.46, -6.42],
        ['C-H', -0.95, 6.15, -6.42],
        ['C-O (ether)', 0.19, 8.19, -4.82],
        ['C-O (ether)', -0.57, 7.7, -6.68],
        ['C-H', 0.17, 5.05, -5.93],
        ['C-H', 1.15, 8.68, -3.98],
        ['C-H', 0.17, 6.96, -4.85],
        ['C-O (ether)', -1.68, 5.62, -7.3],
        ['C-C', -0.62, 6.21, -5.47],
        ['C-H', -0.36, 4.65, -5.68],
        ['C-H', 1.62, 6.55, -3.22],
        ['C-H', 0.12, 5.9, -4.34],
        ['C-C', 0.19, 8.12, -5.43],
        ['C-H', 1.77, 7.5, -4.81],
        ['C-O (alcohol)', 0.81, 8.35, -6.89],
        ['C-C', 0.98, 9.13, -4.56],
        ['C-H', 1.67, 8.54, -3.83],
        ['C-H', -0.81, 5.51, -5.69],
        ['C-C', -1.19, 6.35, -7.24],
        ['C-H', -0.97, 6.06, -6.18],
        ['C-O (alcohol)', -0.59, 7.24, -5.49],
        ['C-C', 0.99, 9.41, -4.99],
        ['C-H', 1.36, 7.0, -2.64],
        ['C-H', 2.62, 9.08, -2.16],
        ['C-C', -0.42, 6.0, -4.79],
        ['C-H', 2.84, 9.42, -5.37],
        ['C-O (alcohol)', -2.93, 5.5, -7.99],
        ['C-H', 0.27, 6.6, -3.87],
        ['C-H', 0.46, 5.44, -4.61],
        ['C-H', 1.47, 6.38, -2.8],
        ['O-H', -0.76, 12.57, -5.81],
        ['O-H', -4.19, -2.06, -5.53],
        ['O-H', -2.67, 2.0, -6.27]
    ])
    sorted_data = custom_sort(data)
    data_bis = sort_Mode([
        ['C-H', 1.18, 8.18, -5.68],
        ['C-H', 2.14, 8.84, -5.27],
        ['C-O (carbonate)', -0.11, 6.69, -6.51],
        ['C-O (carbonate)', -2.54, 7.41, -9.6],
        ['C-H', 0.22, 5.36, -4.08],
        ['C-H', 0.56, 5.2, -3.81],
        ['C-H', 1.18, 5.64, -3.76],
        ['C-O (carbonate)', -1.11, 4.1, -7.36],
        ['C-H', 0.31, 5.54, -3.24],
        ['C-H', 0.7, 6.02, -3.92],
        ['C-H', 1.17, 4.61, -3.49],
        ['C-O (carbonate)', -1.97, 5.84, -8.22],
        ['C-O (carbonate)', -2.04, 4.49, -6.48],
        ['C-O (carbonate)', -2.64, 4.22, -7.51],
        ['C=O (carbonate)', 0.63, 4.24, -2.15],
        ['C-O (carbonate)', -1.6, 5.63, -7.16],
        ['C-O (carbonate)', -2.06, 3.26, -6.69],
        ['C=O (carbonate)', 0.32, 4.94, -2.7],
        ['C-H', 2.59, 5.58, -0.99],
        ['C-H', 2.01, 5.5, -1.14],
        ['C-H', -0.33, 4.0, -2.93],
        ['C-O (urethane)', -1.46, 6.1, -6.81],
        ['C-H', 1.21, 6.42, -3.32],
        ['C-H', 1.31, 7.12, -2.82],
        ['C-N (urethane)', -0.56, 7.1, -5.78],
        ['C-O (urethane)', -0.97, 12.98, -6.8],
        ['C-N (urethane)', 0.97, 8.24, -4.4],
        ['C-O (urethane)', -0.73, 7.63, -6.65],
        ['C=O (urethane)', -0.81, 1.93, -3.62],
        ['C-N (urethane)', -0.03, 4.55, -5.08],
        ['C-O (urethane)', -3.27, 3.39, -7.27],
        ['C=O (urethane)', 0.33, 5.32, -2.8],
        ['C-H', -1.78, 3.84, -6.86],
        ['C-H', 1.08, 6.75, -4.94],
        ['C-H', -0.84, 4.81, -6.16],
        ['C-N (urethane)', -0.27, 5.74, -5.3],
        ['N-H (urethane)', -0.07, 6.09, -4.62],
        ['N-H (urethane)', 1.34, 4.99, -2.87],
        ['C-H', 2.33, 8.2, -3.42],
        ['C-H', 0.72, 6.15, -4.61],
        ['C-H', 0.29, 7.29, -5.1],
        ['C-N (amine)', 0.59, 7.61, -4.6],
        ['C-H', -0.35, 8.58, -6.69],
        ['C-H', -2.27, 6.07, -8.7],
        ['C-N (amine)', -0.23, 6.21, -5.39],
        ['C-N (amine)', 0.61, 9.24, -6.3],
        ['C-H', -0.64, 6.66, -6.74],
        ['C-H', 1.24, 8.64, -4.94],
        ['C-N (amine)', 1.39, 6.33, -5.87],
        ['C-N (amine)', 0.73, 7.55, -5.7],
        ['C-H', -1.42, 3.38, -5.76],
        ['C-H', 0.75, 4.93, -3.7],
        ['C-H', -2.49, 2.39, -6.74],
        ['C-N (amine)', 1.66, 8.66, -3.95],
        ['N-H (amine)', 3.78, 8.02, -2.19],
        ['N-H (amine)', -1.7, 3.73, -3.76],
        ['N-H (amine)', 0.74, 4.87, -2.57],
        ['C-H', 0.72, 6.33, -4.78],
        ['C-H', 0.06, 5.05, -4.27],
        ['C-H', 2.99, 7.76, -1.86],
        ['C-O (ester)', -0.55, 5.9, -7.66],
        ['C-C', -0.58, 7.32, -6.67],
        ['C-H', 2.39, 6.97, -3.82],
        ['C-H', 1.75, 6.56, -3.16],
        ['C-O (ester)', -0.84, 5.66, -7.46],
        ['C-O (ester)', -3.52, 5.24, -8.55],
        ['C=O (ester)', 1.09, 4.39, -1.89],
        ['C-C', 1.27, 12.53, -4.65],
        ['C-H', 0.43, 7.94, -5.85],
        ['C-H', 0.09, 7.74, -5.81],
        ['C-O (ester)', -0.39, 11.23, -7.66],
        ['C-O (ester)', -2.66, 19.02, -8.59],
        ['C=O (ester)', 0.31, 11.62, -2.26],
        ['C-C', 0.24, 6.45, -6.63],
        ['C-H', -0.09, 5.57, -5.87],
        ['C-H', -0.88, 5.48, -5.4],
        ['C-H', 2.44, 7.48, -2.71],
        ['C-O (ester)', -3.88, 2.22, -8.57],
        ['C=O (ester)', 1.46, 5.54, -1.82],
        ['C-H', 1.88, 6.92, -3.07],
        ['C-H', -2.01, 2.66, -6.24],
        ['C-H', -1.17, 4.5, -6.51],
        ['C-O (ether)', -0.99, 5.47, -5.77],
        ['C-H', -0.08, 5.03, -6.98],
        ['C-H', -0.9, 5.57, -6.08],
        ['C-O (ether)', -2.43, 5.22, -8.08],
        ['C-O (ether)', 2.23, 11.05, -6.41],
        ['C-H', -0.36, 5.13, -5.96],
        ['C-H', -0.56, 5.08, -4.56],
        ['C-O (ether)', -0.42, 7.2, -7.31],
        ['C-O (ether)', 1.25, 9.5, -6.22],
        ['C-H', -0.05, 5.22, -4.62],
        ['C-H', 0.53, 8.13, -5.56],
        ['C-H', 0.16, 6.86, -5.42],
        ['C-O (ether)', -0.09, 8.63, -6.02],
        ['C-C', -0.57, 8.46, -6.98],
        ['C-H', 0.25, 10.68, -8.12],
        ['C-H', 0.11, 9.64, -6.93],
        ['C-H', -0.2, 9.37, -6.72],
        ['C-C', 0.41, 9.81, -4.99],
        ['C-H', 1.82, 11.79, -6.17],
        ['C-O (alcohol)', -0.36, 8.82, -8.67],
        ['C-C', 0.82, 9.43, -6.08],
        ['C-H', 0.76, 8.84, -6.26],
        ['C-H', 1.07, 9.26, -6.96],
        ['C-C', -1.08, 6.96, -8.35],
        ['C-H', 0.97, 8.85, -5.1],
        ['C-O (alcohol)', -1.6, 7.87, -7.79],
        ['C-C', 0.34, 11.72, -6.2],
        ['C-H', 2.04, 8.53, -4.58],
        ['C-H', 0.48, 8.07, -5.14],
        ['C-C', -0.22, 7.55, -6.97],
        ['C-H', 0.92, 8.35, -4.35],
        ['C-O (alcohol)', -1.11, 9.38, -9.18],
        ['C-H', 0.63, 7.91, -5.7],
        ['C-H', 0.96, 6.67, -5.35],
        ['C-H', 0.46, 7.3, -5.26],
        ['O-H', -0.38, 100.94, -7.43],
        ['O-H', -2.94, -0.72, -4.17],
        ['O-H', -0.18, 45.39, -5.53]
    ])
    sorted_data_bis = custom_sort(data_bis)
    third_data = sort_Mode([
        ['C-H', 1.14, 8.37, -4.42],
        ['C-H', 0.45, 6.75, -4.4],
        ['C-O (carbonate)', -0.94, 6.13, -7.08],
        ['C-O (carbonate)', -1.4, 7.0, -6.88],
        ['C-H', -0.23, 3.94, -4.16],
        ['C-H', 2.97, 7.28, -1.61],
        ['C-H', 4.37, 8.45, -0.94],
        ['C-O (carbonate)', -1.45, 7.1, -6.8],
        ['C-H', 0.23, 5.97, -4.96],
        ['C-H', -0.3, 5.19, -4.35],
        ['C-H', -0.85, 4.14, -5.37],
        ['C-O (carbonate)', -1.27, 4.92, -7.58],
        ['C-O (carbonate)', -3.61, 3.04, -7.95],
        ['C-O (carbonate)', -1.05, 4.96, -6.95],
        ['C=O (carbonate)', 1.26, 4.29, -1.59],
        ['C-O (carbonate)', -4.56, 0.64, -9.08],
        ['C-O (carbonate)', -1.3, 4.41, -5.88],
        ['C=O (carbonate)', 0.44, 3.13, -2.78],
        ['C-H', 3.64, 10.81, -2.1],
        ['C-H', 0.15, 6.59, -5.68],
        ['C-H', -1.7, 5.16, -7.65],
        ['C-O (urethane)', -1.98, 5.76, -8.7],
        ['C-H', 0.22, 6.24, -4.22],
        ['C-H', -0.28, 5.79, -4.78],
        ['C-N (urethane)', 0.48, 9.56, -4.29],
        ['C-O (urethane)', -0.48, 10.86, -6.2],
        ['C-N (urethane)', -7.67, 0.0, -12.26],
        ['C-O (urethane)', 4.11, 14.12, -5.69],
        ['C=O (urethane)', 3.11, 8.51, -3.11],
        ['C-N (urethane)', -0.76, 3.98, -5.41],
        ['C-O (urethane)', -1.94, 4.06, -8.3],
        ['C=O (urethane)', 1.71, 5.03, -2.42],
        ['C-H', -0.2, 4.12, -6.3],
        ['C-H', -1.39, 3.78, -5.53],
        ['C-H', 1.85, 6.8, -3.62],
        ['C-N (urethane)', 0.14, 7.12, -5.69],
        ['N-H (urethane)', 2.89, 347.15, -2.36],
        ['N-H (urethane)', -1.27, 3.11, -5.07],
        ['C-H', 0.56, 7.52, -5.49],
        ['C-H', -1.17, 4.73, -6.72],
        ['C-H', 2.74, 11.52, -3.76],
        ['C-N (amine)', -0.89, 7.34, -6.97],
        ['C-H', -2.62, 6.03, -9.84],
        ['C-H', -0.35, 8.52, -7.46],
        ['C-N (amine)', -0.08, 9.42, -8.34],
        ['C-N (amine)', 0.97, 9.15, -6.74],
        ['C-H', 0.85, 9.73, -7.85],
        ['C-H', 0.92, 9.38, -6.14],
        ['C-N (amine)', 0.22, 8.2, -6.45],
        ['C-N (amine)', -0.2, 7.49, -6.61],
        ['C-H', -0.19, 4.49, -5.8],
        ['C-H', -0.91, 3.24, -4.81],
        ['C-H', -0.89, 4.47, -6.15],
        ['C-N (amine)', 0.85, 6.9, -5.33],
        ['N-H (amine)', -1.35, 5.38, -7.95],
        ['N-H (amine)', 0.95, 4.45, -3.78],
        ['N-H (amine)', 1.07, 5.35, -3.29],
        ['C-H', 0.15, 4.68, -3.96],
        ['C-H', -0.01, 4.23, -4.34],
        ['C-H', 1.53, 5.43, -2.29],
        ['C-O (ester)', -1.34, 5.9, -6.48],
        ['C-C', 0.49, 8.17, -6.18],
        ['C-H', 1.31, 8.78, -4.0],
        ['C-H', 0.82, 7.69, -3.3],
        ['C-O (ester)', -2.09, 5.61, -8.01],
        ['C-O (ester)', -2.45, 4.06, -7.45],
        ['C=O (ester)', 0.24, 4.2, -2.45],
        ['C-C', -0.51, 6.75, -5.88],
        ['C-H', 2.36, 7.64, -2.82],
        ['C-H', 1.48, 10.04, -4.11],
        ['C-O (ester)', -2.4, 4.05, -8.24],
        ['C-O (ester)', -0.55, 9.07, -7.1],
        ['C=O (ester)', 0.45, 3.52, -3.22],
        ['C-C', -1.66, 5.22, -6.54],
        ['C-H', -1.71, 4.13, -6.59],
        ['C-H', 1.91, 6.05, -3.05],
        ['C-H', 1.54, 7.11, -4.15],
        ['C-O (ester)', 0.68, 13.59, -5.83],
        ['C=O (ester)', 0.04, 2.53, -3.06],
        ['C-H', -2.17, 2.79, -7.45],
        ['C-H', -1.47, 4.82, -5.57],
        ['C-H', 1.93, 9.71, -2.94],
        ['C-O (ether)', -0.12, 8.76, -6.42],
        ['C-H', -1.0, 6.25, -6.36],
        ['C-H', -0.61, 6.9, -5.73],
        ['C-O (ether)', -2.87, 8.29, -9.57],
        ['C-O (ether)', -0.15, 9.76, -6.62],
        ['C-H', -2.95, 3.57, -9.34],
        ['C-H', 1.63, 8.64, -6.12],
        ['C-O (ether)', -1.22, 9.3, -7.01],
        ['C-O (ether)', -2.0, 7.84, -8.42],
        ['C-H', -1.46, 5.0, -6.43],
        ['C-H', 0.07, 5.31, -5.5],
        ['C-H', 1.95, 8.02, -3.11],
        ['C-O (ether)', 0.28, 8.76, -5.03],
        ['C-C', -0.43, 6.56, -7.28],
        ['C-H', 0.0, 7.21, -5.79],
        ['C-H', 0.55, 6.4, -3.95],
        ['C-H', 0.49, 7.5, -4.69],
        ['C-C', 0.24, 8.03, -5.75],
        ['C-H', 0.87, 9.48, -5.5],
        ['C-O (alcohol)', 0.34, 9.8, -5.84],
        ['C-C', 0.39, 9.89, -4.84],
        ['C-H', 1.23, 8.76, -6.14],
        ['C-H', 1.37, 7.94, -5.6],
        ['C-C', -1.3, 8.19, -7.5],
        ['C-H', 1.3, 7.93, -5.2],
        ['C-O (alcohol)', -1.8, 5.51, -8.42],
        ['C-C', 0.05, 7.32, -5.55],
        ['C-H', 1.21, 8.04, -5.75],
        ['C-H', -0.13, 8.2, -5.34],
        ['C-C', -0.49, 7.16, -6.88],
        ['C-H', 0.95, 11.06, -6.52],
        ['C-O (alcohol)', 0.77, 12.46, -6.86],
        ['C-H', 0.52, 6.41, -5.03],
        ['C-H', 0.68, 6.64, -4.94],
        ['C-H', 1.08, 7.95, -4.66],
        ['O-H', -1.23, 8.14, -3.5],
        ['O-H', -1.8, 6.8, -6.33],
        ['O-H', -0.82, 4.37, -4.49]
    ])
    sorted_third_data = custom_sort(third_data)
    fourth_data = sort_Mode([
        ['C-H', 0.73, 7.12, -5.06],
        ['C-H', 0.13, 6.52, -4.79],
        ['C-O (carbonate)', -1.94, 7.35, -8.81],
        ['C-O (carbonate)', -1.34, 11.04, -7.3],
        ['C-H', -1.43, 2.65, -5.77],
        ['C-H', 1.41, 6.39, -2.43],
        ['C-H', -2.34, 1.38, -5.99],
        ['C-O (carbonate)', -1.74, 4.98, -7.83],
        ['C-H', 3.03, 6.67, -0.4],
        ['C-H', -1.7, 1.82, -5.45],
        ['C-H', 1.24, 4.52, -1.97],
        ['C-O (carbonate)', -1.78, 4.26, -7.47],
        ['C-O (carbonate)', -1.54, 4.91, -5.93],
        ['C-O (carbonate)', -0.36, 5.37, -5.97],
        ['C=O (carbonate)', -0.78, 2.78, -3.88],
        ['C-O (carbonate)', -1.31, 5.64, -8.13],
        ['C-O (carbonate)', -1.34, 5.57, -6.45],
        ['C=O (carbonate)', 0.07, 3.08, -2.88],
        ['C-H', 0.27, 6.08, -3.36],
        ['C-H', 0.5, 5.35, -3.56],
        ['C-H', 0.57, 5.59, -3.29],
        ['C-O (urethane)', -1.69, 5.71, -6.93],
        ['C-H', 0.81, 7.84, -4.6],
        ['C-H', 0.56, 6.31, -4.64],
        ['C-N (urethane)', 0.0, 7.11, -5.05],
        ['C-O (urethane)', -1.65, 7.55, -9.66],
        ['C-N (urethane)', 0.91, 7.03, -3.43],
        ['C-O (urethane)', -1.07, 5.51, -5.37],
        ['C=O (urethane)', -0.34, 2.25, -3.21],
        ['C-N (urethane)', -0.12, 5.43, -4.23],
        ['C-O (urethane)', -1.1, 6.75, -6.45],
        ['C=O (urethane)', -0.51, 3.99, -3.32],
        ['C-H', 0.51, 5.63, -4.32],
        ['C-H', 0.66, 6.5, -4.21],
        ['C-H', 0.33, 6.08, -4.28],
        ['C-N (urethane)', -0.44, 4.81, -5.59],
        ['N-H (urethane)', 1.05, 4.27, -3.49],
        ['N-H (urethane)', 0.22, 4.53, -2.48],
        ['C-H', 2.06, 9.34, -3.78],
        ['C-H', -0.7, 5.19, -6.48],
        ['C-H', -1.18, 7.08, -7.03],
        ['C-N (amine)', -0.11, 8.06, -5.66],
        ['C-H', -2.19, 2.85, -7.76],
        ['C-H', -1.82, 5.33, -8.03],
        ['C-N (amine)', -1.71, 7.32, -8.38],
        ['C-N (amine)', 0.11, 9.91, -6.84],
        ['C-H', 0.6, 7.54, -4.32],
        ['C-H', -0.11, 4.27, -4.77],
        ['C-N (amine)', -0.76, 8.83, -6.28],
        ['C-N (amine)', -1.13, 5.02, -8.8],
        ['C-H', 0.84, 7.82, -7.19],
        ['C-H', 1.21, 8.56, -4.22],
        ['C-H', 0.5, 7.19, -5.32],
        ['C-N (amine)', -0.22, 5.87, -5.1],
        ['N-H (amine)', 0.98, 5.82, -3.58],
        ['N-H (amine)', -1.11, 4.0, -5.21],
        ['N-H (amine)', -2.55, 1.59, -6.3],
        ['C-H', -0.57, 3.53, -3.99],
        ['C-H', -0.23, 4.36, -3.1],
        ['C-H', 2.38, 6.46, -0.88],
        ['C-O (ester)', -2.31, 4.48, -7.98],
        ['C-C', -0.34, 7.14, -5.43],
        ['C-H', 0.25, 5.73, -3.49],
        ['C-H', 0.0, 4.91, -3.93],
        ['C-O (ester)', -0.62, 4.33, -6.79],
        ['C-O (ester)', -2.2, 3.96, -6.05],
        ['C=O (ester)', -0.07, 3.33, -2.38],
        ['C-C', 0.66, 6.9, -5.45],
        ['C-H', 1.92, 8.81, -3.38],
        ['C-H', -0.6, 4.41, -6.36],
        ['C-O (ester)', -1.71, 3.71, -7.25],
        ['C-O (ester)', -1.36, 5.88, -7.4],
        ['C=O (ester)', 0.21, 3.46, -2.2],
        ['C-C', -0.01, 6.17, -5.38],
        ['C-H', -0.3, 4.62, -3.94],
        ['C-H', 1.07, 5.19, -1.98],
        ['C-H', -1.32, 1.9, -4.22],
        ['C-O (ester)', -0.84, 7.31, -5.94],
        ['C=O (ester)', -0.61, 2.94, -4.22],
        ['C-H', -2.07, 1.81, -4.82],
        ['C-H', -0.33, 3.9, -3.72],
        ['C-H', 0.45, 3.8, -1.74],
        ['C-O (ether)', -0.98, 4.62, -6.26],
        ['C-H', 1.45, 7.7, -3.57],
        ['C-H', -1.52, 3.77, -5.94],
        ['C-O (ether)', -2.94, 4.17, -8.88],
        ['C-O (ether)', 0.31, 8.89, -5.06],
        ['C-H', 1.5, 5.51, -1.46],
        ['C-H', 0.13, 3.17, -3.64],
        ['C-O (ether)', -1.81, 6.05, -5.95],
        ['C-O (ether)', -0.69, 4.55, -6.26],
        ['C-H', -1.58, 1.12, -4.44],
        ['C-H', 0.56, 2.94, -2.8],
        ['C-H', -0.23, 3.19, -4.16],
        ['C-O (ether)', -1.01, 2.57, -5.31],
        ['C-C', -0.51, 6.08, -5.19],
        ['C-H', -1.16, 3.2, -5.11],
        ['C-H', -0.92, 4.12, -4.49],
        ['C-H', 0.81, 5.12, -3.09],
        ['C-C', 0.19, 8.36, -5.6],
        ['C-H', 1.12, 7.82, -6.0],
        ['C-O (alcohol)', -0.57, 7.97, -6.39],
        ['C-C', 0.26, 7.8, -6.3],
        ['C-H', 0.08, 6.69, -5.9],
        ['C-H', 1.24, 8.04, -4.66],
        ['C-C', -0.14, 7.78, -6.24],
        ['C-H', 0.8, 6.57, -4.87],
        ['C-O (alcohol)', -1.4, 4.35, -6.82],
        ['C-C', 0.61, 8.12, -5.41],
        ['C-H', 2.64, 8.57, -3.42],
        ['C-H', -0.54, 6.62, -5.41],
        ['C-C', 1.69, 7.33, -5.35],
        ['C-H', 1.29, 7.25, -6.46],
        ['C-O (alcohol)', -1.67, 7.44, -8.68],
        ['C-H', -0.88, 4.03, -5.5],
        ['C-H', 1.62, 7.85, -2.99],
        ['C-H', 0.05, 5.69, -4.47],
        ['O-H', 0.21, 6.67, -2.89],
        ['O-H', -2.31, -0.91, -3.8],
        ['O-H', -2.6, -0.11, -4.09]
    ])
    sorted_fourth_data = custom_sort(fourth_data)

    # EF - old version
    # NMC811_EF = [
    # ['C-H', 1.15, 4.88, -3.9] ,
    # ['C-H', 1.09, 5.99, -2.84] ,
    # ['C-O (carbonate)', -1.89, 6.26, -7.91] ,
    # ['C-O (carbonate)', -2.41, 6.06, -7.28] ,
    # ['C-H', 0.48, 5.3, -3.39] ,
    # ['C-H', 0.66, 4.48, -2.95] ,
    # ['C-H', 0.63, 5.05, -3.06] ,
    # ['C-O (carbonate)', -2.58, 3.33, -7.12] ,
    # ['C-H', 0.68, 4.68, -3.73] ,
    # ['C-H', 0.47, 4.48, -3.78] ,
    # ['C-H', 0.39, 5.53, -2.96] ,
    # ['C-O (carbonate)', -1.93, 3.84, -9.14] ,
    # ['C-O (carbonate)', -0.93, 4.99, -5.49] ,
    # ['C-O (carbonate)', -0.82, 6.28, -6.51] ,
    # ['C=O (carbonate)', -0.54, 2.59, -4.14] ,
    # ['C-O (carbonate)', -0.46, 5.33, -6.37] ,
    # ['C-O (carbonate)', -0.63, 3.79, -5.14] ,
    # ['C=O (carbonate)', -0.56, 2.11, -2.98] ,
    # ['C-H', 0.8, 4.12, -2.88] ,
    # ['C-H', 0.1, 4.95, -3.07] ,
    # ['C-H', 0.56, 5.22, -3.33] ,
    # ['C-O (urethane)', -1.55, 3.8, -7.0] ,
    # ['C-H', 0.72, 6.87, -3.64] ,
    # ['C-H', 0.13, 6.72, -3.99] ,
    # ['C-N (urethane)', 0.02, 7.59, -4.72] ,
    # ['C-O (urethane)', -1.89, 5.9, -8.14] ,
    # ['C-N (urethane)', 0.15, 6.47, -4.14] ,
    # ['C-O (urethane)', -0.53, 4.51, -6.1] ,
    # ['C=O (urethane)', -0.65, 2.75, -3.72] ,
    # ['C-N (urethane)', 0.86, 5.24, -4.18] ,
    # ['C-O (urethane)', -1.54, 7.7, -5.04] ,
    # ['C=O (urethane)', -0.77, 2.32, -4.33] ,
    # ['C-H', 0.09, 5.11, -4.09] ,
    # ['C-H', 1.03, 5.85, -3.89] ,
    # ['C-H', 0.41, 5.38, -4.09] ,
    # ['C-N (urethane)', -0.75, 5.5, -5.37] ,
    # ['N-H (urethane)', 0.33, 3.51, -2.91] ,
    # ['N-H (urethane)', 0.22, 5.7, -3.37] ,
    # ['C-H', 0.77, 8.58, -5.25] ,
    # ['C-H', 0.94, 7.91, -5.12] ,
    # ['C-H', -0.58, 11.63, -5.68] ,
    # ['C-N (amine)', -0.96, 5.08, -7.65] ,
    # ['C-H', 0.25, 5.52, -5.35] ,
    # ['C-H', 0.28, 6.05, -5.16] ,
    # ['C-N (amine)', -1.13, 7.43, -7.2] ,
    # ['C-N (amine)', -1.14, 8.1, -6.22] ,
    # ['C-H', 0.71, 9.06, -4.75] ,
    # ['C-H', 0.31, 7.7, -5.61] ,
    # ['C-N (amine)', -1.98, 6.87, -6.85] ,
    # ['C-N (amine)', -2.03, 5.51, -7.16] ,
    # ['C-H', 0.28, 10.32, -5.25] ,
    # ['C-H', 0.22, 6.2, -5.19] ,
    # ['C-H', 1.22, 6.82, -4.97] ,
    # ['C-N (amine)', -1.73, 4.69, -8.54] ,
    # ['N-H (amine)', -0.18, 5.63, -4.31] ,
    # ['N-H (amine)', -0.41, 4.73, -3.24] ,
    # ['N-H (amine)', -0.65, 5.88, -5.48] ,
    # ['C-H', -0.1, 3.35, -2.45] ,
    # ['C-H', 0.98, 4.54, -2.36] ,
    # ['C-H', 0.28, 4.0, -2.59] ,
    # ['C-O (ester)', -2.14, 3.04, -7.42] ,
    # ['C-C', -0.04, 4.75, -4.54] ,
    # ['C-H', 0.67, 4.23, -2.55] ,
    # ['C-H', 0.62, 3.68, -2.63] ,
    # ['C-O (ester)', -1.64, 5.16, -6.26] ,
    # ['C-O (ester)', -2.13, 2.01, -5.87] ,
    # ['C=O (ester)', 0.37, 3.03, -1.51] ,
    # ['C-C', 0.64, 8.17, -5.15] ,
    # ['C-H', 0.81, 5.14, -3.79] ,
    # ['C-H', 0.53, 5.95, -3.82] ,
    # ['C-O (ester)', -2.11, 4.71, -7.36] ,
    # ['C-O (ester)', -0.88, 6.64, -5.2] ,
    # ['C=O (ester)', -0.22, 3.06, -2.75] ,
    # ['C-C', -0.34, 4.26, -4.84] ,
    # ['C-H', 0.68, 3.61, -2.83] ,
    # ['C-H', 0.75, 3.7, -1.62] ,
    # ['C-H', 0.3, 3.3, -2.02] ,
    # ['C-O (ester)', -0.17, 5.96, -7.81] ,
    # ['C=O (ester)', -0.31, 2.41, -2.58] ,
    # ['C-H', 1.14, 5.33, -2.36] ,
    # ['C-H', 0.08, 4.08, -3.13] ,
    # ['C-H', 0.69, 4.99, -3.36] ,
    # ['C-O (ether)', -1.14, 4.19, -6.79] ,
    # ['C-H', -0.02, 4.01, -3.97] ,
    # ['C-H', -0.06, 5.38, -4.43] ,
    # ['C-O (ether)', -1.44, 5.08, -6.07] ,
    # ['C-O (ether)', -1.54, 8.55, -7.35] ,
    # ['C-H', -0.02, 4.56, -5.97] ,
    # ['C-H', -0.19, 6.08, -5.02] ,
    # ['C-O (ether)', -1.03, 8.46, -6.25] ,
    # ['C-O (ether)', -0.55, 5.58, -6.81] ,
    # ['C-H', 0.45, 3.88, -4.19] ,
    # ['C-H', 1.12, 4.96, -2.49] ,
    # ['C-H', 0.15, 4.48, -2.93] ,
    # ['C-O (ether)', -1.47, 6.81, -6.87] ,
    # ['C-C', -0.75, 7.51, -5.91] ,
    # ['C-H', 0.28, 5.65, -4.4] ,
    # ['C-H', 0.14, 7.03, -4.95] ,
    # ['C-H', 0.85, 5.12, -4.07] ,
    # ['C-C', -0.21, 8.92, -5.24] ,
    # ['C-H', 1.41, 8.0, -5.59] ,
    # ['C-O (alcohol)', -0.75, 10.71, -7.33] ,
    # ['C-C', 0.7, 8.22, -5.84] ,
    # ['C-H', 1.19, 7.47, -5.22] ,
    # ['C-H', -0.1, 8.17, -4.85] ,
    # ['C-C', -0.74, 5.4, -7.26] ,
    # ['C-H', 2.11, 7.68, -4.57] ,
    # ['C-O (alcohol)', -1.52, 5.0, -6.97] ,
    # ['C-C', 1.2, 7.09, -5.05] ,
    # ['C-H', 0.74, 8.0, -5.08] ,
    # ['C-H', 0.56, 7.44, -5.66] ,
    # ['C-C', 0.27, 6.32, -5.67] ,
    # ['C-H', 0.81, 7.22, -3.73] ,
    # ['C-O (alcohol)', -2.26, 3.97, -7.47] ,
    # ['C-H', 0.25, 5.15, -3.37] ,
    # ['C-H', 0.39, 3.96, -3.1] ,
    # ['C-H', 0.76, 7.1, -4.31] ,
    # ['O-H', -1.52, 9.14, -5.77] ,
    # ['O-H', -2.97, 0.44, -4.69] ,
    # ['O-H', -2.1, 1.33, -3.98]
    # ]
    NMC811_EF = [['C-H', 0.8, 4.29, -3.05],
                 ['C-H', 1.37, 5.56, -2.42],
                 ['C-O (carbonate)', -0.36, 10.36, -6.52],
                 ['C-O (carbonate)', -1.88, 7.71, -8.15],
                 ['C-H', 0.2, 5.93, -3.86],
                 ['C-H', 0.24, 5.77, -3.84],
                 ['C-H', 0.64, 5.52, -3.77],
                 ['C-O (carbonate)', -1.5, 4.08, -8.73],
                 ['C-H', 0.61, 4.58, -2.45],
                 ['C-H', 0.5, 4.34, -2.71],
                 ['C-H', 0.75, 3.46, -2.38],
                 ['C-O (carbonate)', -1.76, 4.37, -5.94],
                 ['C-O (carbonate)', -0.63, 5.89, -5.86],
                 ['C-O (carbonate)', -2.21, 5.01, -6.74],
                 ['C=O (carbonate)', -0.59, 2.96, -3.01],
                 ['C-O (carbonate)', -0.62, 6.79, -5.88],
                 ['C-O (carbonate)', -0.94, 3.4, -5.94],
                 ['C=O (carbonate)', -0.74, 1.93, -2.96],
                 ['C-H', 0.46, 4.47, -3.32],
                 ['C-H', 0.69, 4.81, -3.2],
                 ['C-H', 0.77, 4.42, -3.43],
                 ['C-O (urethane)', -1.55, 5.03, -6.83],
                 ['C-H', 1.24, 7.52, -4.67],
                 ['C-H', 1.11, 8.56, -4.63],
                 ['C-N (urethane)', 0.15, 6.33, -5.67],
                 ['C-O (urethane)', -2.92, 4.87, -7.39],
                 ['C-N (urethane)', 0.63, 7.2, -4.64],
                 ['C-O (urethane)', -1.98, 4.01, -6.06],
                 ['C=O (urethane)', -0.22, 2.75, -2.51],
                 ['C-N (urethane)', 0.17, 5.2, -4.97],
                 ['C-O (urethane)', 0.3, 7.72, -7.58],
                 ['C=O (urethane)', -0.78, 3.06, -3.4],
                 ['C-H', 0.09, 5.36, -4.68],
                 ['C-H', 0.75, 5.78, -4.18],
                 ['C-H', -0.17, 7.35, -4.87],
                 ['C-N (urethane)', -0.58, 4.64, -5.49],
                 ['N-H (urethane)', 0.71, 4.24, -2.62],
                 ['N-H (urethane)', 0.93, 4.31, -2.53],
                 ['C-H', 0.82, 5.96, -4.63],
                 ['C-H', 0.85, 7.84, -5.22],
                 ['C-H', -0.83, 6.85, -5.24],
                 ['C-N (amine)', -0.45, 5.81, -5.47],
                 ['C-H', 1.54, 8.21, -4.74],
                 ['C-H', 0.51, 7.55, -5.86],
                 ['C-N (amine)', -1.84, 5.63, -6.92],
                 ['C-N (amine)', 1.69, 10.27, -5.18],
                 ['C-H', 1.43, 7.33, -4.12],
                 ['C-H', 0.57, 6.79, -5.05],
                 ['C-N (amine)', -1.72, 4.71, -6.48],
                 ['C-N (amine)', -0.13, 5.92, -5.39],
                 ['C-H', 1.15, 5.66, -5.22],
                 ['C-H', 0.65, 5.65, -3.36],
                 ['C-H', 0.06, 6.98, -3.42],
                 ['C-N (amine)', -1.08, 3.88, -6.55],
                 ['N-H (amine)', -0.5, 2.86, -2.99],
                 ['N-H (amine)', 0.06, 4.45, -3.75],
                 ['N-H (amine)', -0.84, 4.55, -3.81],
                 ['C-H', 0.18, 3.89, -2.46],
                 ['C-H', 0.75, 4.25, -2.42],
                 ['C-H', 0.58, 3.57, -2.58],
                 ['C-O (ester)', -1.58, 5.83, -6.84],
                 ['C-C', 0.55, 8.03, -4.75],
                 ['C-H', 0.66, 5.17, -4.34],
                 ['C-H', 0.78, 5.93, -3.52],
                 ['C-O (ester)', -1.29, 4.88, -5.5],
                 ['C-O (ester)', -2.82, 2.59, -6.66],
                 ['C=O (ester)', 0.2, 3.12, -2.78],
                 ['C-C', 1.33, 7.8, -5.83],
                 ['C-H', 0.3, 5.41, -3.18],
                 ['C-H', 1.02, 5.18, -3.57],
                 ['C-O (ester)', -1.09, 5.2, -7.02],
                 ['C-O (ester)', -2.61, 2.79, -7.9],
                 ['C=O (ester)', 0.35, 3.96, -2.55],
                 ['C-C', 0.08, 4.2, -4.85],
                 ['C-H', 0.64, 5.3, -3.93],
                 ['C-H', 0.78, 5.26, -2.81],
                 ['C-H', 0.23, 4.85, -3.33],
                 ['C-O (ester)', -1.45, 4.54, -6.26],
                 ['C=O (ester)', 0.25, 3.04, -1.99],
                 ['C-H', 0.81, 7.17, -4.08],
                 ['C-H', 0.13, 5.79, -4.6],
                 ['C-H', -0.1, 4.67, -3.93],
                 ['C-O (ether)', -0.73, 4.48, -6.72],
                 ['C-H', -0.62, 5.4, -5.26],
                 ['C-H', 0.03, 5.84, -5.03],
                 ['C-O (ether)', -1.73, 5.98, -7.2],
                 ['C-O (ether)', -0.39, 8.66, -7.08],
                 ['C-H', 0.23, 5.74, -4.96],
                 ['C-H', -0.62, 5.73, -5.32],
                 ['C-O (ether)', -1.52, 11.94, -6.62],
                 ['C-O (ether)', -1.23, 5.59, -7.4],
                 ['C-H', -0.12, 4.68, -4.12],
                 ['C-H', 0.86, 6.32, -3.51],
                 ['C-H', -0.17, 5.04, -4.69],
                 ['C-O (ether)', -1.19, 3.91, -5.62],
                 ['C-C', 0.0, 7.84, -6.28],
                 ['C-H', 0.64, 6.1, -3.73],
                 ['C-H', -0.21, 5.17, -4.18],
                 ['C-H', 1.11, 6.45, -4.56],
                 ['C-C', 0.17, 7.28, -4.67],
                 ['C-H', 1.64, 7.91, -4.26],
                 ['C-O (alcohol)', -0.59, 9.3, -7.91],
                 ['C-C', 0.01, 7.52, -4.69],
                 ['C-H', 1.29, 8.66, -4.38],
                 ['C-H', 1.12, 6.19, -4.27],
                 ['C-C', -1.25, 6.82, -8.54],
                 ['C-H', 3.23, 8.27, -3.86],
                 ['C-O (alcohol)', -1.79, 8.27, -7.57],
                 ['C-C', 0.11, 8.19, -4.94],
                 ['C-H', 1.41, 8.17, -4.02],
                 ['C-H', 0.67, 5.66, -4.04],
                 ['C-C', 0.04, 6.5, -5.57],
                 ['C-H', 1.62, 8.28, -3.31],
                 ['C-O (alcohol)', -1.69, 4.67, -7.46],
                 ['C-H', 0.39, 5.41, -3.98],
                 ['C-H', 0.48, 5.86, -4.21],
                 ['C-H', 0.49, 5.35, -3.4],
                 ['O-H', -1.7, 3.79, -4.0],
                 ['O-H', -2.37, 2.42, -4.56],
                 ['O-H', -2.21, 0.3, -3.37]]
    sorted_NMC811_EF = custom_sort(NMC811_EF)

    # EF - old version
    # NMC622_EF = [
    # ['C-H', 0.75, 4.6, -2.88] ,
    # ['C-H', 1.01, 6.01, -3.15] ,
    # ['C-O (carbonate)', -1.09, 6.28, -6.11] ,
    # ['C-O (carbonate)', -1.88, 5.16, -7.75] ,
    # ['C-H', 0.63, 4.23, -2.89] ,
    # ['C-H', 0.53, 4.35, -3.37] ,
    # ['C-H', 0.33, 4.45, -3.11] ,
    # ['C-O (carbonate)', -2.27, 3.14, -7.16] ,
    # ['C-H', 0.48, 4.02, -2.6] ,
    # ['C-H', 0.52, 4.45, -2.63] ,
    # ['C-H', 0.66, 3.96, -2.21] ,
    # ['C-O (carbonate)', -1.68, 2.95, -7.23] ,
    # ['C-O (carbonate)', -0.84, 4.13, -5.34] ,
    # ['C-O (carbonate)', -0.62, 3.27, -5.88] ,
    # ['C=O (carbonate)', -0.53, 1.86, -2.69] ,
    # ['C-O (carbonate)', -1.26, 3.88, -6.1] ,
    # ['C-O (carbonate)', -0.72, 5.01, -4.61] ,
    # ['C=O (carbonate)', -0.58, 2.29, -2.98] ,
    # ['C-H', 0.43, 4.51, -2.71] ,
    # ['C-H', 0.61, 4.32, -3.03] ,
    # ['C-H', 0.59, 4.38, -2.83] ,
    # ['C-O (urethane)', -2.1, 3.1, -5.94] ,
    # ['C-H', 0.54, 5.56, -2.87] ,
    # ['C-H', 0.71, 6.24, -3.38] ,
    # ['C-N (urethane)', -0.45, 5.16, -5.51] ,
    # ['C-O (urethane)', -1.72, 7.98, -7.54] ,
    # ['C-N (urethane)', 0.48, 5.99, -3.77] ,
    # ['C-O (urethane)', -0.53, 5.14, -5.66] ,
    # ['C=O (urethane)', -0.74, 2.58, -3.63] ,
    # ['C-N (urethane)', 0.4, 6.29, -3.45] ,
    # ['C-O (urethane)', -0.62, 6.7, -5.4] ,
    # ['C=O (urethane)', -0.89, 2.12, -3.84] ,
    # ['C-H', 0.32, 4.92, -3.84] ,
    # ['C-H', 0.86, 5.52, -4.28] ,
    # ['C-H', 0.12, 6.12, -4.32] ,
    # ['C-N (urethane)', -0.88, 4.35, -5.44] ,
    # ['N-H (urethane)', -0.11, 3.9, -2.72] ,
    # ['N-H (urethane)', 0.49, 4.26, -3.2] ,
    # ['C-H', 0.84, 8.08, -4.51] ,
    # ['C-H', -0.77, 4.16, -5.05] ,
    # ['C-H', 0.48, 7.86, -4.34] ,
    # ['C-N (amine)', -1.26, 5.67, -5.79] ,
    # ['C-H', 0.33, 7.25, -6.45] ,
    # ['C-H', 0.33, 7.5, -6.15] ,
    # ['C-N (amine)', -1.7, 7.48, -7.2] ,
    # ['C-N (amine)', -0.9, 7.57, -6.43] ,
    # ['C-H', 1.16, 9.13, -4.98] ,
    # ['C-H', 1.17, 7.13, -5.67] ,
    # ['C-N (amine)', -0.28, 7.52, -5.81] ,
    # ['C-N (amine)', -1.2, 5.25, -7.52] ,
    # ['C-H', -0.38, 4.96, -6.05] ,
    # ['C-H', 0.7, 7.11, -4.65] ,
    # ['C-H', 0.14, 7.33, -4.75] ,
    # ['C-N (amine)', -0.85, 6.47, -6.27] ,
    # ['N-H (amine)', -0.06, 4.98, -4.43] ,
    # ['N-H (amine)', 0.43, 5.91, -5.63] ,
    # ['N-H (amine)', 0.03, 3.72, -4.95] ,
    # ['C-H', 0.2, 6.17, -4.03] ,
    # ['C-H', 0.95, 7.21, -4.49] ,
    # ['C-H', 1.05, 6.36, -4.16] ,
    # ['C-O (ester)', -2.0, 5.23, -7.8] ,
    # ['C-C', 1.49, 8.11, -5.01] ,
    # ['C-H', 0.65, 5.44, -4.04] ,
    # ['C-H', 1.35, 6.4, -4.23] ,
    # ['C-O (ester)', -2.37, 6.48, -7.69] ,
    # ['C-O (ester)', -1.53, 6.13, -7.2] ,
    # ['C=O (ester)', -0.2, 2.96, -2.63] ,
    # ['C-C', 0.32, 10.27, -4.47] ,
    # ['C-H', 0.75, 7.31, -4.96] ,
    # ['C-H', 0.45, 6.74, -4.59] ,
    # ['C-O (ester)', -1.59, 5.56, -8.06] ,
    # ['C-O (ester)', -2.16, 5.32, -7.05] ,
    # ['C=O (ester)', 0.2, 3.86, -3.54] ,
    # ['C-C', -0.74, 6.24, -6.56] ,
    # ['C-H', 0.6, 6.21, -5.13] ,
    # ['C-H', 0.37, 6.86, -4.43] ,
    # ['C-H', 0.66, 5.51, -3.72] ,
    # ['C-O (ester)', -1.27, 5.49, -6.91] ,
    # ['C=O (ester)', -0.53, 3.83, -2.92] ,
    # ['C-H', -0.09, 4.84, -4.2] ,
    # ['C-H', -0.22, 5.32, -3.69] ,
    # ['C-H', 1.19, 6.14, -2.94] ,
    # ['C-O (ether)', -1.19, 2.81, -5.43] ,
    # ['C-H', -0.64, 3.78, -4.19] ,
    # ['C-H', -0.27, 3.42, -3.74] ,
    # ['C-O (ether)', -1.09, 3.91, -5.79] ,
    # ['C-O (ether)', -0.59, 5.53, -5.49] ,
    # ['C-H', -0.48, 2.82, -3.85] ,
    # ['C-H', -0.38, 2.9, -3.46] ,
    # ['C-O (ether)', -0.22, 4.97, -5.52] ,
    # ['C-O (ether)', -1.21, 4.17, -5.9] ,
    # ['C-H', 0.5, 4.66, -3.93] ,
    # ['C-H', 0.51, 5.66, -3.62] ,
    # ['C-H', 0.13, 4.44, -3.63] ,
    # ['C-O (ether)', -0.81, 3.09, -5.29] ,
    # ['C-C', -0.84, 6.28, -5.69] ,
    # ['C-H', 0.92, 5.92, -3.98] ,
    # ['C-H', 0.42, 4.89, -4.53] ,
    # ['C-H', -0.06, 6.21, -3.82] ,
    # ['C-C', 0.14, 8.16, -4.88] ,
    # ['C-H', 1.02, 8.62, -5.32] ,
    # ['C-O (alcohol)', -0.75, 9.68, -7.28] ,
    # ['C-C', 0.87, 8.52, -4.55] ,
    # ['C-H', 0.83, 7.44, -4.26] ,
    # ['C-H', -0.13, 7.45, -5.07] ,
    # ['C-C', -1.09, 6.87, -7.74] ,
    # ['C-H', 1.35, 9.14, -5.15] ,
    # ['C-O (alcohol)', -0.84, 5.59, -6.67] ,
    # ['C-C', -0.19, 7.07, -5.56] ,
    # ['C-H', 0.9, 6.81, -4.66] ,
    # ['C-H', 1.26, 7.49, -4.38] ,
    # ['C-C', -0.01, 7.33, -5.01] ,
    # ['C-H', 0.76, 6.86, -4.53] ,
    # ['C-O (alcohol)', -2.72, 5.42, -8.26] ,
    # ['C-H', 0.09, 5.85, -4.53] ,
    # ['C-H', 0.16, 5.07, -3.95] ,
    # ['C-H', 0.88, 5.16, -3.86] ,
    # ['O-H', -1.12, 2.4, -3.73] ,
    # ['O-H', -2.82, 1.17, -5.32] ,
    # ['O-H', -2.15, 0.01, -3.5]
    # ]
    NMC622_EF = [['C-H', 0.83, 6.28, -4.74],
                 ['C-H', 0.91, 7.96, -3.9],
                 ['C-O (carbonate)', 0.82, 7.07, -7.99],
                 ['C-O (carbonate)', -2.91, 4.77, -9.87],
                 ['C-H', 0.79, 5.02, -3.69],
                 ['C-H', 0.43, 5.13, -3.33],
                 ['C-H', 0.54, 4.58, -2.7],
                 ['C-O (carbonate)', -1.51, 3.98, -8.09],
                 ['C-H', 0.37, 3.84, -2.19],
                 ['C-H', 0.33, 3.46, -2.17],
                 ['C-H', 0.37, 3.74, -2.93],
                 ['C-O (carbonate)', -1.54, 3.76, -6.54],
                 ['C-O (carbonate)', -0.81, 4.21, -5.18],
                 ['C-O (carbonate)', -1.33, 5.86, -7.1],
                 ['C=O (carbonate)', -0.48, 3.11, -2.54],
                 ['C-O (carbonate)', -0.35, 8.32, -6.75],
                 ['C-O (carbonate)', -1.37, 4.32, -5.93],
                 ['C=O (carbonate)', -0.39, 2.67, -2.96],
                 ['C-H', 0.95, 5.14, -3.31],
                 ['C-H', 1.01, 4.85, -3.07],
                 ['C-H', 0.71, 5.14, -2.98],
                 ['C-O (urethane)', -3.03, 3.83, -9.49],
                 ['C-H', 1.25, 5.91, -3.13],
                 ['C-H', 0.66, 5.62, -3.65],
                 ['C-N (urethane)', 0.96, 8.57, -4.67],
                 ['C-O (urethane)', -2.91, 5.93, -10.0],
                 ['C-N (urethane)', -0.25, 6.26, -5.26],
                 ['C-O (urethane)', -0.2, 6.75, -7.01],
                 ['C=O (urethane)', -0.25, 3.51, -4.24],
                 ['C-N (urethane)', -0.55, 5.32, -4.54],
                 ['C-O (urethane)', 0.07, 7.58, -6.19],
                 ['C=O (urethane)', -0.45, 3.62, -3.22],
                 ['C-H', 0.49, 4.84, -3.78],
                 ['C-H', 1.08, 5.49, -3.43],
                 ['C-H', 0.28, 4.72, -3.9],
                 ['C-N (urethane)', -0.71, 5.83, -5.95],
                 ['N-H (urethane)', 0.38, 7.23, -4.09],
                 ['N-H (urethane)', 0.64, 4.54, -2.75],
                 ['C-H', 0.11, 5.27, -4.03],
                 ['C-H', 0.37, 5.63, -4.59],
                 ['C-H', 0.62, 6.51, -3.07],
                 ['C-N (amine)', -0.89, 4.89, -5.16],
                 ['C-H', -0.52, 6.31, -4.8],
                 ['C-H', 0.28, 4.63, -4.34],
                 ['C-N (amine)', -0.63, 6.81, -6.36],
                 ['C-N (amine)', -0.46, 6.08, -5.91],
                 ['C-H', 1.73, 7.05, -3.48],
                 ['C-H', 0.74, 5.26, -4.47],
                 ['C-N (amine)', -1.42, 4.81, -5.84],
                 ['C-N (amine)', -1.12, 4.69, -5.07],
                 ['C-H', 0.32, 6.08, -5.1],
                 ['C-H', 0.3, 8.91, -4.03],
                 ['C-H', 0.64, 6.65, -3.81],
                 ['C-N (amine)', -0.5, 5.0, -5.29],
                 ['N-H (amine)', 0.29, 3.63, -3.89],
                 ['N-H (amine)', -0.21, 4.89, -3.97],
                 ['N-H (amine)', -0.22, 2.92, -2.85],
                 ['C-H', 0.2, 3.98, -3.08],
                 ['C-H', 0.28, 4.13, -3.35],
                 ['C-H', 0.81, 5.27, -2.22],
                 ['C-O (ester)', -0.98, 6.21, -6.74],
                 ['C-C', 0.75, 9.35, -4.86],
                 ['C-H', 0.85, 5.55, -3.77],
                 ['C-H', 0.27, 5.39, -3.51],
                 ['C=O (ester)', -1.54, 5.39, -7.33],
                 ['C-O (ester)', -2.55, 3.22, -7.78],
                 ['C-O (ester)', 0.7, 4.01, -1.74],
                 ['C-C', 0.94, 7.09, -4.2],
                 ['C-H', 0.69, 5.26, -2.55],
                 ['C-H', 0.41, 3.9, -2.87],
                 ['C-O (ester)', -1.74, 5.0, -6.93],
                 ['C-O (ester)', -0.6, 7.16, -6.34],
                 ['C=O (ester)', -0.03, 3.04, -2.12],
                 ['C-C', -0.37, 4.72, -4.98],
                 ['C-H', 0.15, 3.3, -2.08],
                 ['C-H', 0.61, 3.61, -1.7],
                 ['C-H', 0.59, 3.6, -2.25],
                 ['C-O (ester)', -1.28, 3.99, -5.87],
                 ['C=O (ester)', -0.34, 2.77, -2.54],
                 ['C-H', 0.04, 4.15, -4.33],
                 ['C-H', 0.09, 5.15, -4.49],
                 ['C-H', 0.92, 6.89, -3.49],
                 ['C-O (ether)', -0.29, 5.56, -5.58],
                 ['C-H', -0.13, 5.72, -6.33],
                 ['C-H', -0.33, 6.72, -5.54],
                 ['C-O (ether)', -2.28, 2.33, -7.18],
                 ['C-O (ether)', 0.17, 8.01, -5.63],
                 ['C-H', -0.1, 5.38, -4.86],
                 ['C-H', -0.78, 4.29, -5.56],
                 ['C-O (ether)', -1.78, 5.13, -7.32],
                 ['C-O (ether)', 0.72, 10.14, -5.18],
                 ['C-H', 0.18, 4.97, -5.11],
                 ['C-H', 1.67, 5.76, -3.44],
                 ['C-H', 1.16, 4.82, -3.58],
                 ['C-O (ether)', -2.29, 3.39, -6.84],
                 ['C-C', -0.67, 6.68, -6.8],
                 ['C-H', 1.03, 6.35, -4.18],
                 ['C-H', 0.66, 5.5, -3.73],
                 ['C-H', 0.32, 5.71, -4.37],
                 ['C-C', 0.51, 8.71, -6.34],
                 ['C-H', 2.06, 8.9, -3.65],
                 ['C-O (alcohol)', -2.27, 4.84, -7.98],
                 ['C-C', 0.68, 7.56, -5.06],
                 ['C-H', 1.2, 9.65, -5.31],
                 ['C-H', 0.82, 10.38, -6.26],
                 ['C-C', -1.76, 6.63, -7.02],
                 ['C-H', 2.24, 6.32, -3.28],
                 ['C-O (alcohol)', 0.21, 6.61, -6.58],
                 ['C-C', 0.36, 8.3, -5.75],
                 ['C-H', 1.31, 8.39, -5.64],
                 ['C-H', 0.14, 7.28, -5.14],
                 ['C-C', -0.13, 6.89, -5.1],
                 ['C-H', 0.66, 7.9, -3.69],
                 ['C-O (alcohol)', -3.13, 5.44, -8.73],
                 ['C-H', 0.62, 5.81, -4.09],
                 ['C-H', 0.27, 6.45, -3.59],
                 ['C-H', 0.68, 5.76, -3.73],
                 ['O-H', -1.6, 3.25, -3.73],
                 ['O-H', -3.18, -0.56, -5.11],
                 ['O-H', -1.64, 7.3, -4.86]]
    sorted_NMC622_EF = custom_sort(NMC622_EF)

    # EF - old version
    # LiNMC811_EF = [
    # ['C-H', 0.96, 7.61, -4.16] ,
    # ['C-H', 1.53, 7.69, -4.5] ,
    # ['C-O (carbonate)', -1.05, 8.86, -8.68] ,
    # ['C-O (carbonate)', -0.89, 9.05, -8.86] ,
    # ['C-H', 0.38, 5.7, -4.43] ,
    # ['C-H', 0.57, 5.23, -4.18] ,
    # ['C-H', 0.33, 6.0, -4.04] ,
    # ['C-O (carbonate)', -2.26, 6.26, -7.37] ,
    # ['C-H', 0.79, 5.93, -4.46] ,
    # ['C-H', 0.51, 6.15, -4.2] ,
    # ['C-H', 0.29, 7.22, -5.05] ,
    # ['C-O (carbonate)', -2.18, 6.35, -7.99] ,
    # ['C-O (carbonate)', -1.58, 3.86, -6.91] ,
    # ['C-O (carbonate)', -2.76, 6.83, -7.59] ,
    # ['C=O (carbonate)', 1.05, 4.3, -1.9] ,
    # ['C-O (carbonate)', -1.79, 4.71, -7.91] ,
    # ['C-O (carbonate)', -1.9, 5.62, -7.23] ,
    # ['C=O (carbonate)', 0.36, 4.18, -2.96] ,
    # ['C-H', 0.91, 4.97, -3.4] ,
    # ['C-H', 0.11, 4.41, -3.47] ,
    # ['C-H', 0.85, 5.22, -3.43] ,
    # ['C-O (urethane)', -1.65, 5.43, -6.46] ,
    # ['C-H', 1.36, 8.73, -5.31] ,
    # ['C-H', 0.27, 8.31, -5.44] ,
    # ['C-N (urethane)', -0.52, 6.68, -6.53] ,
    # ['C-O (urethane)', -1.47, 7.82, -7.79] ,
    # ['C-N (urethane)', 0.29, 6.51, -4.12] ,
    # ['C-O (urethane)', -0.42, 7.37, -5.94] ,
    # ['C=O (urethane)', -0.85, 1.92, -3.67] ,
    # ['C-N (urethane)', -0.02, 4.11, -4.22] ,
    # ['C-O (urethane)', -2.25, 4.68, -7.17] ,
    # ['C=O (urethane)', 0.17, 3.7, -2.55] ,
    # ['C-H', 0.33, 4.44, -3.4] ,
    # ['C-H', 0.22, 4.34, -2.71] ,
    # ['C-H', 0.21, 4.29, -2.8] ,
    # ['C-N (urethane)', -0.63, 5.84, -5.21] ,
    # ['N-H (urethane)', 1.61, 14.9, -3.43] ,
    # ['N-H (urethane)', 0.15, 4.91, -3.82] ,
    # ['C-H', 1.18, 6.37, -3.9] ,
    # ['C-H', 0.4, 9.22, -4.45] ,
    # ['C-H', -0.28, 6.38, -6.06] ,
    # ['C-N (amine)', -0.91, 5.68, -5.79] ,
    # ['C-H', 0.23, 8.46, -5.85] ,
    # ['C-H', 0.87, 7.51, -5.13] ,
    # ['C-N (amine)', 0.03, 6.04, -6.35] ,
    # ['C-N (amine)', -0.92, 7.72, -6.11] ,
    # ['C-H', 0.39, 8.72, -5.05] ,
    # ['C-H', -0.12, 6.71, -5.37] ,
    # ['C-N (amine)', -1.76, 5.11, -7.97] ,
    # ['C-N (amine)', 1.69, 10.31, -4.69] ,
    # ['C-H', -0.2, 5.36, -6.1] ,
    # ['C-H', 1.28, 7.57, -4.19] ,
    # ['C-H', 0.52, 5.8, -4.45] ,
    # ['C-N (amine)', 0.28, 7.09, -6.14] ,
    # ['N-H (amine)', -0.01, 5.26, -4.06] ,
    # ['N-H (amine)', 0.23, 4.12, -3.28] ,
    # ['N-H (amine)', -0.35, 5.5, -3.72] ,
    # ['C-H', 0.83, 5.78, -3.64] ,
    # ['C-H', 1.18, 5.8, -3.54] ,
    # ['C-H', 0.45, 4.8, -3.64] ,
    # ['C-O (ester)', -0.56, 8.36, -7.42] ,
    # ['C-C', 0.38, 7.11, -5.25] ,
    # ['C-H', 0.07, 8.01, -4.5] ,
    # ['C-H', 0.2, 8.16, -5.0] ,
    # ['C-O (ester)', -1.22, 4.37, -7.4] ,
    # ['C-O (ester)', -4.04, 3.87, -8.4] ,
    # ['C=O (ester)', 0.96, 5.43, -2.48] ,
    # ['C-C', -0.17, 8.63, -7.81] ,
    # ['C-H', 0.9, 9.67, -6.42] ,
    # ['C-H', 1.71, 12.69, -5.76] ,
    # ['C-O (ester)', -2.49, 7.01, -7.73] ,
    # ['C-O (ester)', -0.58, 10.56, -6.99] ,
    # ['C=O (ester)', -0.16, 4.48, -3.73] ,
    # ['C-C', -0.79, 6.78, -5.93] ,
    # ['C-H', 1.1, 5.61, -4.32] ,
    # ['C-H', 1.25, 7.49, -3.75] ,
    # ['C-H', 0.72, 5.44, -3.72] ,
    # ['C-O (ester)', -1.0, 7.71, -7.44] ,
    # ['C=O (ester)', -0.41, 3.97, -3.64] ,
    # ['C-H', 1.12, 5.64, -3.31] ,
    # ['C-H', 0.57, 4.81, -3.92] ,
    # ['C-H', -0.07, 4.61, -5.04] ,
    # ['C-O (ether)', -0.44, 6.68, -5.48] ,
    # ['C-H', 0.04, 4.97, -4.78] ,
    # ['C-H', -0.56, 4.37, -5.01] ,
    # ['C-O (ether)', -0.05, 9.83, -5.97] ,
    # ['C-O (ether)', -0.85, 6.3, -8.44] ,
    # ['C-H', -0.28, 4.37, -5.21] ,
    # ['C-H', -0.37, 4.11, -5.75] ,
    # ['C-O (ether)', -0.6, 8.03, -6.73] ,
    # ['C-O (ether)', -1.09, 7.95, -5.85] ,
    # ['C-H', 0.39, 5.56, -4.68] ,
    # ['C-H', 0.75, 7.91, -4.98] ,
    # ['C-H', -0.06, 5.66, -4.12] ,
    # ['C-O (ether)', -0.63, 5.81, -5.4] ,
    # ['C-C', -0.51, 7.02, -6.07] ,
    # ['C-H', 0.42, 6.83, -4.94] ,
    # ['C-H', 0.52, 7.87, -5.87] ,
    # ['C-H', 1.05, 5.9, -4.09] ,
    # ['C-C', 0.73, 7.64, -5.03] ,
    # ['C-H', 1.23, 9.12, -5.09] ,
    # ['C-O (alcohol)', -1.15, 6.84, -7.06] ,
    # ['C-C', -0.36, 9.27, -5.42] ,
    # ['C-H', 0.26, 7.7, -5.25] ,
    # ['C-H', 1.01, 8.55, -5.55] ,
    # ['C-C', -1.27, 6.89, -7.06] ,
    # ['C-H', 2.01, 9.94, -4.87] ,
    # ['C-O (alcohol)', 0.02, 8.77, -7.29] ,
    # ['C-C', 0.52, 9.3, -6.24] ,
    # ['C-H', 0.74, 8.39, -5.64] ,
    # ['C-H', 0.68, 8.18, -5.33] ,
    # ['C-C', -0.24, 9.06, -7.11] ,
    # ['C-H', 1.4, 7.6, -4.69] ,
    # ['C-O (alcohol)', -3.05, 7.24, -10.05] ,
    # ['C-H', 0.45, 6.89, -5.14] ,
    # ['C-H', 0.04, 6.22, -5.55] ,
    # ['C-H', 0.11, 6.72, -5.19] ,
    # ['O-H', 0.06, 8.08, -5.4] ,
    # ['O-H', -1.3, 6.28, -5.84] ,
    # ['O-H', -1.12, 7.2, -4.94]
    # ]
    LiNMC811_EF = [['C-H', 1.11, 6.47, -4.35],
                   ['C-H', 0.98, 5.51, -3.65],
                   ['C-O (carbonate)', -0.8, 8.05, -7.9],
                   ['C-O (carbonate)', -2.03, 6.66, -8.2],
                   ['C-H', 0.18, 4.04, -3.18],
                   ['C-H', 0.37, 5.18, -3.48],
                   ['C-H', 0.44, 5.69, -3.94],
                   ['C-O (carbonate)', -0.98, 5.96, -8.97],
                   ['C-H', 0.99, 5.29, -3.76],
                   ['C-H', 0.51, 5.68, -4.16],
                   ['C-H', 0.37, 6.53, -4.01],
                   ['C-O (carbonate)', -0.86, 6.61, -7.09],
                   ['C-O (carbonate)', -1.35, 3.08, -6.24],
                   ['C-O (carbonate)', -1.44, 5.05, -6.85],
                   ['C=O (carbonate)', -0.17, 3.51, -2.1],
                   ['C-O (carbonate)', -1.13, 5.01, -5.86],
                   ['C-O (carbonate)', -2.78, 4.5, -6.12],
                   ['C=O (carbonate)', 0.27, 4.34, -2.73],
                   ['C-H', 1.03, 7.06, -4.64],
                   ['C-H', 0.83, 6.09, -4.54],
                   ['C-H', 1.25, 7.19, -4.58],
                   ['C-O (urethane)', -1.17, 5.3, -7.52],
                   ['C-H', 1.58, 9.71, -5.49],
                   ['C-H', 1.44, 9.82, -6.78],
                   ['C-N (urethane)', -0.04, 8.91, -6.13],
                   ['C-O (urethane)', -2.22, 8.91, -9.19],
                   ['C-N (urethane)', 0.01, 6.67, -4.72],
                   ['C-O (urethane)', -3.39, 2.63, -8.42],
                   ['C=O (urethane)', 1.55, 5.32, -2.97],
                   ['C-N (urethane)', -1.07, 4.86, -5.1],
                   ['C-O (urethane)', -1.83, 7.4, -8.1],
                   ['C=O (urethane)', 0.7, 4.45, -3.61],
                   ['C-H', 0.07, 6.61, -5.61],
                   ['C-H', 0.65, 7.11, -4.96],
                   ['C-H', 0.13, 6.96, -6.05],
                   ['C-N (urethane)', -0.14, 6.89, -5.77],
                   ['N-H (urethane)', 0.59, 5.46, -3.64],
                   ['N-H (urethane)', -0.32, 4.22, -3.82],
                   ['C-H', 0.76, 7.08, -4.51],
                   ['C-H', 0.62, 6.94, -4.24],
                   ['C-H', -0.79, 6.19, -4.97],
                   ['C-N (amine)', 0.29, 5.97, -4.25],
                   ['C-H', -1.0, 6.46, -4.76],
                   ['C-H', -0.79, 7.11, -6.63],
                   ['C-N (amine)', 0.61, 9.33, -5.83],
                   ['C-N (amine)', -0.83, 6.99, -5.94],
                   ['C-H', 0.82, 7.33, -4.84],
                   ['C-H', -1.02, 6.34, -5.43],
                   ['C-N (amine)', -1.74, 6.64, -7.2],
                   ['C-N (amine)', 1.52, 11.53, -3.76],
                   ['C-H', -0.26, 5.29, -5.07],
                   ['C-H', 0.62, 6.93, -3.5],
                   ['C-H', -0.03, 7.1, -5.3],
                   ['C-N (amine)', 0.58, 6.88, -4.47],
                   ['N-H (amine)', -0.06, 6.64, -3.6],
                   ['N-H (amine)', -2.09, 5.1, -4.82],
                   ['N-H (amine)', -0.32, 3.5, -3.39],
                   ['C-H', 0.34, 4.67, -3.58],
                   ['C-H', 0.8, 5.9, -4.2],
                   ['C-H', -0.04, 4.21, -3.6],
                   ['C-O (ester)', -0.6, 5.22, -6.55],
                   ['C-C', 0.78, 6.65, -5.44],
                   ['C-H', 0.64, 5.96, -4.17],
                   ['C-H', 1.08, 6.69, -4.64],
                   ['C-O (ester)', 0.14, 8.02, -6.19],
                   ['C-O (ester)', -3.91, 1.76, -9.25],
                   ['C=O (ester)', 1.26, 5.0, -1.99],
                   ['C-C', 0.23, 7.56, -5.88],
                   ['C-H', 0.34, 8.83, -5.93],
                   ['C-H', 0.65, 7.3, -4.67],
                   ['C-O (ester)', -0.67, 6.43, -6.13],
                   ['C-O (ester)', -3.7, 3.34, -8.1],
                   ['C=O (ester)', 1.45, 4.56, -1.87],
                   ['C-C', -1.32, 6.61, -5.88],
                   ['C-H', -0.12, 5.41, -4.05],
                   ['C-H', 0.7, 7.89, -3.62],
                   ['C-H', 0.08, 5.6, -3.53],
                   ['C-O (ester)', -1.84, 3.97, -8.26],
                   ['C=O (ester)', 0.92, 4.07, -2.3],
                   ['C-H', 0.85, 6.08, -4.27],
                   ['C-H', -0.02, 6.0, -5.04],
                   ['C-H', 0.58, 4.85, -4.53],
                   ['C-O (ether)', -0.59, 6.4, -5.16],
                   ['C-H', 0.02, 5.31, -5.38],
                   ['C-H', -0.18, 5.6, -5.21],
                   ['C-O (ether)', -0.56, 8.0, -7.14],
                   ['C-O (ether)', -1.81, 6.55, -6.21],
                   ['C-H', 0.02, 4.8, -5.26],
                   ['C-H', -0.61, 5.89, -5.22],
                   ['C-O (ether)', -0.82, 10.06, -7.35],
                   ['C-O (ether)', 0.53, 6.79, -6.38],
                   ['C-H', 0.05, 4.12, -4.1],
                   ['C-H', 0.24, 3.57, -2.44],
                   ['C-H', 0.15, 3.92, -3.7],
                   ['C-O (ether)', 0.51, 5.73, -4.61],
                   ['C-C', -0.8, 6.31, -6.5],
                   ['C-H', 0.38, 6.69, -3.73],
                   ['C-H', 0.19, 8.58, -5.29],
                   ['C-H', 0.55, 8.45, -4.6],
                   ['C-C', 0.42, 7.73, -5.45],
                   ['C-H', 1.39, 7.38, -4.03],
                   ['C-O (alcohol)', -0.89, 9.59, -7.18],
                   ['C-C', 0.32, 8.29, -6.08],
                   ['C-H', 0.73, 7.56, -4.97],
                   ['C-H', 0.94, 7.26, -4.77],
                   ['C-C', -0.56, 6.19, -7.17],
                   ['C-H', 1.49, 8.07, -3.19],
                   ['C-O (alcohol)', -1.89, 5.87, -7.63],
                   ['C-C', 0.53, 8.8, -6.1],
                   ['C-H', -0.16, 8.21, -5.39],
                   ['C-H', 0.32, 7.43, -3.95],
                   ['C-C', -0.13, 7.47, -6.47],
                   ['C-H', 0.99, 11.22, -7.22],
                   ['C-O (alcohol)', -0.71, 6.7, -8.27],
                   ['C-H', 0.07, 5.82, -5.13],
                   ['C-H', -0.28, 7.3, -5.02],
                   ['C-H', 0.78, 6.31, -4.41],
                   ['O-H', -0.9, 29.1, -6.83],
                   ['O-H', -2.76, -1.03, -4.82],
                   ['O-H', -1.19, 7.03, -6.34]]
    sorted_LiNMC811_EF = custom_sort(LiNMC811_EF)

    # EF - old version
    # LiNMC622_EF = [
    # ['C-H', 0.8, 8.2, -4.9] ,
    # ['C-H', 1.34, 8.8, -4.89] ,
    # ['C-O (carbonate)', -1.92, 7.49, -8.87] ,
    # ['C-O (carbonate)', -0.68, 7.65, -7.36] ,
    # ['C-H', 0.38, 4.47, -3.04] ,
    # ['C-H', 0.07, 4.61, -2.9] ,
    # ['C-H', 0.24, 4.35, -3.05] ,
    # ['C-O (carbonate)', -1.6, 5.64, -7.16] ,
    # ['C-H', 0.42, 7.66, -4.94] ,
    # ['C-H', 0.94, 7.48, -4.75] ,
    # ['C-H', 0.77, 6.95, -4.82] ,
    # ['C-O (carbonate)', -1.19, 4.27, -6.75] ,
    # ['C-O (carbonate)', -1.66, 4.97, -6.77] ,
    # ['C-O (carbonate)', -2.28, 6.51, -7.34] ,
    # ['C=O (carbonate)', 0.12, 3.96, -2.96] ,
    # ['C-O (carbonate)', -2.54, 4.14, -7.9] ,
    # ['C-O (carbonate)', -1.88, 4.4, -6.72] ,
    # ['C=O (carbonate)', 0.85, 5.16, -4.73] ,
    # ['C-H', 0.74, 5.19, -3.44] ,
    # ['C-H', 0.44, 5.64, -3.1] ,
    # ['C-H', 0.17, 4.65, -3.58] ,
    # ['C-O (urethane)', -1.96, 4.58, -8.11] ,
    # ['C-H', 1.09, 8.04, -4.48] ,
    # ['C-H', 0.61, 7.81, -5.17] ,
    # ['C-N (urethane)', -0.08, 6.76, -7.65] ,
    # ['C-O (urethane)', 0.03, 17.06, -7.44] ,
    # ['C-N (urethane)', 0.57, 6.26, -4.58] ,
    # ['C-O (urethane)', -1.68, 6.69, -6.7] ,
    # ['C=O (urethane)', -0.49, 5.27, -3.52] ,
    # ['C-N (urethane)', -0.4, 5.77, -4.6] ,
    # ['C-O (urethane)', -3.45, 4.91, -8.26] ,
    # ['C=O (urethane)', 0.61, 4.55, -3.1] ,
    # ['C-H', 0.65, 7.0, -5.33] ,
    # ['C-H', 0.49, 6.25, -4.52] ,
    # ['C-H', 0.49, 5.55, -4.47] ,
    # ['C-N (urethane)', 0.21, 6.54, -6.24] ,
    # ['N-H (urethane)', 3.33, 17.19, -3.65] ,
    # ['N-H (urethane)', 0.52, 2.48, -2.24] ,
    # ['C-H', 0.64, 6.16, -3.97] ,
    # ['C-H', 0.26, 4.87, -5.49] ,
    # ['C-H', 1.2, 5.56, -3.43] ,
    # ['C-N (amine)', -1.32, 6.17, -6.62] ,
    # ['C-H', 0.08, 6.25, -4.3] ,
    # ['C-H', 0.82, 6.7, -4.29] ,
    # ['C-N (amine)', -0.95, 7.31, -5.5] ,
    # ['C-N (amine)', 0.02, 6.55, -5.65] ,
    # ['C-H', 0.44, 7.6, -5.3] ,
    # ['C-H', 0.36, 5.98, -4.93] ,
    # ['C-N (amine)', -0.18, 7.65, -6.46] ,
    # ['C-N (amine)', -0.21, 8.57, -6.19] ,
    # ['C-H', -0.61, 4.92, -6.02] ,
    # ['C-H', 0.49, 7.33, -4.99] ,
    # ['C-H', -0.32, 5.56, -4.16] ,
    # ['C-N (amine)', -0.42, 7.22, -6.55] ,
    # ['N-H (amine)', -0.04, 4.36, -3.72] ,
    # ['N-H (amine)', 0.49, 3.63, -4.12] ,
    # ['N-H (amine)', -0.58, 5.61, -3.83] ,
    # ['C-H', 0.42, 5.76, -4.23] ,
    # ['C-H', 0.47, 5.79, -4.54] ,
    # ['C-H', 0.95, 6.34, -3.36] ,
    # ['C-O (ester)', -2.29, 5.5, -6.22] ,
    # ['C-C', -0.12, 7.2, -5.75] ,
    # ['C-H', 0.4, 6.45, -4.44] ,
    # ['C-H', 0.99, 9.44, -3.86] ,
    # ['C-O (ester)', -1.04, 7.25, -6.46] ,
    # ['C-O (ester)', -1.25, 7.29, -5.95] ,
    # ['C=O (ester)', -0.06, 3.31, -3.5] ,
    # ['C-C', 0.38, 6.69, -5.1] ,
    # ['C-H', 0.22, 6.35, -5.11] ,
    # ['C-H', 1.36, 7.36, -4.3] ,
    # ['C-O (ester)', -1.57, 4.21, -6.59] ,
    # ['C-O (ester)', -2.76, 2.25, -8.07] ,
    # ['C=O (ester)', 0.56, 4.06, -2.25] ,
    # ['C-C', -1.28, 5.65, -6.83] ,
    # ['C-H', 0.19, 6.47, -5.05] ,
    # ['C-H', 1.84, 7.88, -2.81] ,
    # ['C-H', 0.39, 7.37, -4.56] ,
    # ['C-O (ester)', -1.66, 4.84, -7.1] ,
    # ['C=O (ester)', 0.99, 4.24, -2.19] ,
    # ['C-H', 0.02, 5.06, -4.54] ,
    # ['C-H', -0.05, 5.56, -4.03] ,
    # ['C-H', 1.06, 5.62, -3.07] ,
    # ['C-O (ether)', -1.27, 4.88, -5.46] ,
    # ['C-H', -0.42, 4.3, -5.24] ,
    # ['C-H', -0.54, 4.47, -3.82] ,
    # ['C-O (ether)', 0.0, 6.39, -6.03] ,
    # ['C-O (ether)', 0.01, 8.36, -6.96] ,
    # ['C-H', -0.16, 5.64, -5.63] ,
    # ['C-H', -0.33, 4.66, -5.03] ,
    # ['C-O (ether)', -1.41, 3.94, -6.87] ,
    # ['C-O (ether)', -0.44, 8.35, -5.28] ,
    # ['C-H', 0.19, 4.84, -4.18] ,
    # ['C-H', 0.72, 6.46, -3.7] ,
    # ['C-H', 0.29, 4.85, -4.74] ,
    # ['C-O (ether)', -0.4, 5.99, -5.95] ,
    # ['C-C', -0.57, 6.97, -6.15] ,
    # ['C-H', 0.85, 7.0, -3.7] ,
    # ['C-H', 0.6, 5.56, -3.89] ,
    # ['C-H', 1.43, 7.05, -4.54] ,
    # ['C-C', -0.26, 8.02, -6.38] ,
    # ['C-H', 1.25, 10.6, -6.88] ,
    # ['C-O (alcohol)', -0.17, 7.64, -7.83] ,
    # ['C-C', 0.7, 8.71, -4.99] ,
    # ['C-H', 0.93, 9.74, -6.13] ,
    # ['C-H', -0.14, 11.44, -7.71] ,
    # ['C-C', -0.79, 7.46, -7.41] ,
    # ['C-H', 0.97, 9.35, -5.75] ,
    # ['C-O (alcohol)', -2.78, 6.18, -8.87] ,
    # ['C-C', -0.42, 7.86, -6.8] ,
    # ['C-H', 0.48, 9.22, -5.47] ,
    # ['C-H', 0.47, 9.31, -5.8] ,
    # ['C-C', 0.18, 6.74, -5.07] ,
    # ['C-H', 3.31, 8.93, -4.29] ,
    # ['C-O (alcohol)', -2.66, 9.07, -10.32] ,
    # ['C-H', 0.21, 6.28, -4.47] ,
    # ['C-H', 1.07, 6.18, -4.35] ,
    # ['C-H', 0.61, 7.5, -4.08] ,
    # ['O-H', 1.13, 7.17, -5.54] ,
    # ['O-H', -3.01, 8.3, -6.29] ,
    # ['O-H', 0.87, 87.67, -6.16]
    # ]
    LiNMC622_EF = [['C-H', 0.9, 8.95, -4.82],
                   ['C-H', 1.52, 8.52, -4.1],
                   ['C-O (carbonate)', -1.39, 6.63, -7.81],
                   ['C-O (carbonate)', -0.95, 8.86, -7.35],
                   ['C-H', 0.54, 5.83, -3.87],
                   ['C-H', 0.58, 5.89, -3.67],
                   ['C-H', 0.4, 5.23, -3.31],
                   ['C-O (carbonate)', 0.32, 6.46, -7.33],
                   ['C-H', 0.44, 8.43, -5.0],
                   ['C-H', 0.98, 7.42, -4.4],
                   ['C-H', 0.84, 7.11, -3.98],
                   ['C-O (carbonate)', -0.98, 6.04, -7.58],
                   ['C-O (carbonate)', -3.37, 2.02, -7.69],
                   ['C-O (carbonate)', -1.33, 6.67, -6.77],
                   ['C=O (carbonate)', 0.83, 4.34, -2.09],
                   ['C-O (carbonate)', -1.44, 5.76, -7.69],
                   ['C-O (carbonate)', -2.33, 5.05, -7.56],
                   ['C=O (carbonate)', 0.7, 4.55, -3.07],
                   ['C-H', 0.23, 7.35, -4.42],
                   ['C-H', 0.88, 5.31, -4.01],
                   ['C-H', 1.22, 6.81, -4.48],
                   ['C-O (urethane)', -1.44, 5.54, -7.5],
                   ['C-H', 1.24, 7.63, -3.69],
                   ['C-H', 1.06, 8.07, -4.12],
                   ['C-N (urethane)', 0.66, 7.97, -5.44],
                   ['C-O (urethane)', -1.83, 6.87, -9.39],
                   ['C-N (urethane)', -1.44, 6.89, -6.66],
                   ['C-O (urethane)', -2.95, 4.39, -7.95],
                   ['C=O (urethane)', 1.08, 7.09, -3.19],
                   ['C-N (urethane)', -0.89, 4.19, -4.64],
                   ['C-O (urethane)', -1.94, 8.25, -7.8],
                   ['C=O (urethane)', 0.78, 4.7, -3.22],
                   ['C-H', 0.58, 6.95, -4.74],
                   ['C-H', 1.49, 7.57, -4.97],
                   ['C-H', 0.96, 6.81, -5.25],
                   ['C-N (urethane)', -0.32, 7.04, -6.72],
                   ['N-H (urethane)', 0.58, 6.29, -3.47],
                   ['N-H (urethane)', 0.15, 4.11, -3.16],
                   ['C-H', 1.09, 7.69, -5.39],
                   ['C-H', 0.83, 7.27, -6.59],
                   ['C-H', 1.38, 10.6, -5.17],
                   ['C-N (amine)', 0.65, 7.01, -5.34],
                   ['C-H', -0.13, 6.14, -4.57],
                   ['C-H', 0.3, 5.97, -4.44],
                   ['C-N (amine)', 0.9, 8.58, -5.26],
                   ['C-N (amine)', -0.08, 5.9, -6.69],
                   ['C-H', 0.1, 6.81, -5.5],
                   ['C-H', 0.7, 6.56, -5.27],
                   ['C-N (amine)', -0.14, 6.71, -6.34],
                   ['C-N (amine)', -0.37, 7.94, -5.57],
                   ['C-H', -0.66, 6.21, -6.17],
                   ['C-H', 0.59, 6.84, -4.84],
                   ['C-H', -0.34, 6.9, -5.06],
                   ['C-N (amine)', 0.35, 7.21, -5.86],
                   ['N-H (amine)', 1.53, 5.05, -2.68],
                   ['N-H (amine)', 0.94, 7.22, -4.09],
                   ['N-H (amine)', 0.45, 5.16, -3.19],
                   ['C-H', 0.56, 5.39, -3.81],
                   ['C-H', 0.71, 5.8, -4.27],
                   ['C-H', 0.52, 5.96, -3.67],
                   ['C-O (ester)', -1.16, 5.91, -6.83],
                   ['C-C', 0.85, 8.14, -4.77],
                   ['C-H', -0.08, 6.39, -5.47],
                   ['C-H', 0.75, 6.89, -5.28],
                   ['C-O (ester)', -1.94, 4.8, -7.35],
                   ['C-O (ester)', -3.25, 2.18, -8.36],
                   ['C=O (ester)', 0.73, 4.78, -2.04],
                   ['C-C', -0.57, 10.22, -6.25],
                   ['C-H', 0.74, 8.35, -4.22],
                   ['C-H', 0.04, 6.47, -4.15],
                   ['C-O (ester)', -0.88, 6.09, -6.47],
                   ['C-O (ester)', -2.37, 5.14, -8.38],
                   ['C=O (ester)', 0.75, 4.28, -2.19],
                   ['C-C', -0.78, 7.21, -6.69],
                   ['C-H', 1.26, 5.94, -5.63],
                   ['C-H', 0.84, 7.35, -3.77],
                   ['C-H', 0.22, 5.3, -4.2],
                   ['C-O (ester)', -2.34, 3.39, -7.97],
                   ['C=O (ester)', 0.06, 4.44, -2.45],
                   ['C-H', 0.23, 4.32, -3.61],
                   ['C-H', -0.29, 4.75, -4.02],
                   ['C-H', 0.87, 4.64, -2.68],
                   ['C-O (ether)', -0.17, 6.18, -5.07],
                   ['C-H', -0.48, 6.97, -5.98],
                   ['C-H', -0.06, 6.57, -6.06],
                   ['C-O (ether)', -0.49, 7.2, -7.17],
                   ['C-O (ether)', -0.94, 12.31, -7.31],
                   ['C-H', -0.48, 4.98, -5.42],
                   ['C-H', 0.06, 5.27, -5.48],
                   ['C-O (ether)', 2.07, 9.83, -5.54],
                   ['C-O (ether)', -1.35, 5.27, -6.89],
                   ['C-H', -0.61, 4.92, -4.87],
                   ['C-H', 1.35, 5.33, -3.4],
                   ['C-H', 0.62, 5.08, -4.73],
                   ['C-O (ether)', -0.27, 5.79, -4.94],
                   ['C-C', -1.14, 5.54, -5.87],
                   ['C-H', 0.93, 6.22, -4.06],
                   ['C-H', 0.33, 5.53, -3.84],
                   ['C-H', 0.26, 5.59, -4.57],
                   ['C-C', 0.69, 6.2, -5.65],
                   ['C-H', 0.57, 7.62, -4.86],
                   ['C-O (alcohol)', -0.13, 7.64, -6.61],
                   ['C-C', 1.02, 8.88, -5.81],
                   ['C-H', 0.54, 7.23, -5.28],
                   ['C-H', 1.43, 9.71, -5.48],
                   ['C-C', -1.39, 6.72, -7.39],
                   ['C-H', 0.22, 8.86, -4.54],
                   ['C-O (alcohol)', -1.48, 9.56, -7.83],
                   ['C-C', 1.14, 8.93, -6.49],
                   ['C-H', 0.38, 8.24, -5.19],
                   ['C-H', 0.36, 8.77, -6.67],
                   ['C-C', 1.29, 8.7, -5.69],
                   ['C-H', 2.81, 9.9, -3.87],
                   ['C-O (alcohol)', -4.35, 4.15, -11.16],
                   ['C-H', 0.3, 6.51, -5.32],
                   ['C-H', -0.11, 5.91, -4.56],
                   ['C-H', 0.39, 6.63, -4.04],
                   ['O-H', -0.26, 9.76, -5.09],
                   ['O-H', -2.93, -0.56, -4.84],
                   ['O-H', 324.0, 439.13, 174.45]
                   ]
    sorted_LiNMC622_EF = custom_sort(LiNMC622_EF)

    names_data = {'LiNMC811': data_bis, 'LiNMC622': third_data,
                  'NMC811': data, 'NMC622': fourth_data,
                  'LiNMC811 EF': LiNMC811_EF, 'LiNMC622 EF': LiNMC622_EF,
                  'NMC811 EF': NMC811_EF, 'NMC622 EF': NMC622_EF}

    comp = input_comp()

    if comp.upper() == "YES":
        comp_type = input_comp_type()

        materials = append_materials(comp, comp_type, efOrNot=None, several_data=None, names_data=names_data)
        material_data_dico = append_mat_dico(materials, names_data)

        draw_comparative_subplots(comp_type, material_data_dico)
        plot_jointplots(material_data_dico, comp_type)
    else:
        efOrNot = input_data()
        several_data = input_several_data()
        materials = append_materials(comp, comp_type= None, efOrNot=efOrNot, several_data=several_data, names_data=names_data)
        material_data_dico = append_mat_dico(materials, names_data)

        if several_data == 1 and len(materials) == 1:
            mat = materials[0]
            ordered_set, sorted_names_list, sorted_min_values_list, sorted_max_values_list, sorted_mode_values_list = sort_data(several_data, materials, names_data)
            draw_subplot_single_mat(mat, material_data_dico, ordered_set, sorted_names_list, sorted_min_values_list,
                                    sorted_max_values_list)
        else:
            draw_several_subplots(materials, names_data)