from tkinter import filedialog, simpledialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, AutoLocator

def process_excel(input_file, sheet_name, ignore_columns_indices, end_column_index):
    # Read the Excel file with specified sheet name and columns
    df = pd.read_excel(input_file, sheet_name=sheet_name, usecols=range(end_column_index))

    # Drop the columns specified in ignore_columns_indices
    df.drop(df.columns[ignore_columns_indices], axis=1, inplace=True)
    df['Series'] = df['Electrode'].astype(str) + ' ' + df['Molecule'].astype(str)

    return df

def draw_distances(df, barWidth):
    molecule_to_dist = {'Alcohol': [], 'Amine': [],'Carbonate': [], 'Ester': [], 'Ether': [], 'Urethane': []}
    molecule_to_chrg = {'Alcohol': [], 'Amine': [],'Carbonate': [], 'Ester': [], 'Ether': [], 'Urethane': []}

    # Use colors from mcolors.TABLEAU_COLORS
    colors = [mcolors.CSS4_COLORS['darkturquoise'], mcolors.CSS4_COLORS['darkviolet'], mcolors.CSS4_COLORS['red'],
              mcolors.CSS4_COLORS['royalblue'], mcolors.CSS4_COLORS['forestgreen'], mcolors.CSS4_COLORS['gold']]

    # Convert the 'Distance' column to numeric, setting non-numeric values to NaN
    df['Distance'] = pd.to_numeric(df['Distance'], errors='coerce')
    grouped = df.groupby(['Electrode', 'Molecule'], sort=False)

    # Iterate through each group
    for (electrode, molecule), group in grouped:
        #print(group)
        # Filter out NaN values and find the minimum distance
        min_distance = group['Distance'].dropna().min()

        if molecule in molecule_to_dist and pd.notna(min_distance):
            molecule_to_dist[molecule].append(min_distance)
            min_index = group['Distance'].dropna().idxmin()
            molecule_to_chrg[molecule].append(group['Charge'][min_index])
            #print('Distance :', min_distance, 'and charge :', group['Charge'][min_index])



    bar1 = np.arange(len(molecule_to_dist['Alcohol']))
    bar2 = [x + barWidth for x in bar1]
    bar3 = [x + barWidth for x in bar2]
    bar4 = [x + barWidth for x in bar3]
    bar5 = [x + barWidth for x in bar4]
    bar6 = [x + barWidth for x in bar5]
    brs = [bar1, bar2, bar3, bar4, bar5, bar6]

    for molecule, bar, col in zip(molecule_to_dist.keys(), brs, colors):
        print(molecule)
        print('Charge :', molecule_to_chrg[molecule])
        print('bar :', bar)
        print('Distance :', molecule_to_dist[molecule])
        plt.bar(bar, molecule_to_dist[molecule], color=col, width=barWidth, edgecolor='grey', label=molecule)

    # Adding Xticks
    plt.title('Distance to electrode')
    plt.xlabel('Electrode material', fontsize=25)#fontweight='bold'
    plt.ylabel('Distance (\u00C5)', fontsize=25)#fontweight='bold'
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.gca().set_ylim(top=5.5)
    plt.gca().yaxis.set_major_locator(MultipleLocator(1))

    # Calculate the midpoint of each x-category and apply the xtick positions
    #xtick_positions = [r + half_bar_group_width for r in range(len(molecule_to_dist['Alcohol']))]
    xtick_positions = [start + barWidth * 2.5 for start in bar1]
    plt.xticks(xtick_positions, ['LiNMC811', 'LiNMC622', 'NMC811', 'NMC622']) # Adapt

    def draw_bader_charge():
        num_categories = len(molecule_to_dist['Alcohol'])
        num_molecules = len(molecule_to_chrg)

        ax2 = plt.gca().twinx()
        scatter_color = mcolors.CSS4_COLORS['chocolate']  # Choose a color for the scatter points
        x_positions = []
        charges = []

        for i, molecule in enumerate(molecule_to_chrg.keys()):
            # Iterate through each category in molecule_to_chrg[molecule]
            for category_index, charge in enumerate(molecule_to_chrg[molecule]):
                # Calculate the x-position using the list of bar positions (`brs`) for the current molecule
                # and category index to get the correct x position
                x_position = brs[i][category_index]

                # Append the calculated x-position and charge
                x_positions.append(x_position)
                charges.append(charge)

        # Plot the scatter points on the secondary y-axis
        ax2.scatter(x_positions, charges, facecolors='white', edgecolors=scatter_color , label='Bader net atomic charge', s=75)
        ax2.spines['right'].set_color(scatter_color)
        ax2.set_ylabel('Bader net atomic charge', color=scatter_color, fontsize=25) #fontweight='bold'
        ax2.tick_params(axis='y', labelcolor=scatter_color, labelsize=25)

        ax2.invert_yaxis()
        ax2.set_ylim(top=-1.35)
        ax2.set_ylim(bottom=-0.95)
        ax2.yaxis.set_major_locator(MultipleLocator(0.1))

        # Display the legend for the scatter plot
        #ax2.legend()

    draw_bader_charge()

    #plt.legend()
    plt.show()



if __name__ == "__main__":
    ignored_column = [2,3] # Adapt (for my file, these were 'Atom' and 'Index'
    last_column = 6 # Adapt to your system -> last column not taken into account and count from 0
    barWidth = 0.15
    # Read data from Excel sheet
    input_file = filedialog.askopenfilename(title="Select the Excel file")
    sheet_name = simpledialog.askstring("Sheet name", "What is the name of the sheet that contains the data to plot?")
    data = pd.read_excel(input_file, sheet_name=sheet_name)
    df = process_excel(input_file, sheet_name, ignored_column, last_column)
    draw_distances(df, barWidth)
