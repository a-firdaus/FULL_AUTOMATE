#Copyright Lieven Bekaert 2021,2022,2023 all rights reserved - lieven.bekaert@vub.be
print("Script Copyright LIEVEN BEKAERT 2021,22,23 All rights reserved - lieven.bekaert@vub.be")
print("\n\nThis script calculates the RBC for any VASP Generated POSCAR/OUTCAR structure")
#Run pip install wheel pandas numpy matplotlib vpython pyautogui
import warnings
import matplotlib
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from vpython import vector as vpythvector
from vpython import dot, degrees, diff_angle
import tkinter as tk
from tkinter import filedialog
from copy import deepcopy
from molatomclasses import atom
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from copy import deepcopy
from pyautogui import confirm
from matplotlib.lines import Line2D
import csv
matplotlib.use('TkAgg')
warnings.filterwarnings("ignore")


radiuslist={"H":1.20,"Li":1.82,"Cl":2.75 ,"C":1.70, "Si":2.10 ,"N":1.55, "O":1.52,"P":1.80,"S":1.80,'Ne':1.54,'Ar':1.88,'He':1.40,"Na":2.27,"F":1.47,"Br":1.85,"Tc":2.05,"Pm":2.43,"Cd":1.55,"Co":1.35,"Mn":1.61,"Fe":1.61,"Ni":1.49}#IN ANGSTROM (Molecular coord also in A) (from http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html)
colourdict = {"H": (0.8, 0.8, 0.8), "Li":(252/255, 144/255, 245/255),"Cl": (0, 1, 0), "C": (0, 0, 0), "O": (1, 0, 0), "N": (0, 0, 1), "Ar": (0.90, 0.90, 0.90),
              "Ne": (0.90, 0.90, 0.90), "S": (1, 1, 0), "He": (0.90, 0.90, 0.90), "Na": (1, 1, 1), "Pm": (1, 0.8, 1),
              "Tc": (0.8, 1, 1), "F": (159 / 255.0, 252 / 255.0, 206 / 255.0), "Br": (153 / 255.0, 34 / 255.0, 0),
              "Cd": (0.8, 1, 1), "Co": (0.8, 1, 1),"Si":(153/255, 51/255, 242/255),"Ni":(0.5, 1, 1),"Mn":(0.3, 1, 1),"P":(180/255.0, 52/255.0, 235/255.0),"Fe":(0.3, 1, 1)}
molarweightlist={"H":1.00794,"O":15.9994,"C":12.0107,"Ar":39.948,"He":4.002602,"Ne":20.1797,"N":14.0067,"Cl":35.453,"Na":22.989769,"F":18.998403,"Br":79.904,"Tc":98.0,"Pm":145.0,"Cd":112.411,"Co":58.9332,"fulvic_acid":308.242,"Pb":207.2,'Cu':63.546,"Fe":55.845}
element_colors = {'H': 'oldlace', 'Li': 'green', 'O': 'red', 'C': 'saddlebrown', 'N': 'blue',
                      'Mn': 'darkorchid', 'Ni': 'silver', 'Co': 'indigo'}

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
def euclidiandistance_vectors(pos1,pos2):
    return sqrt(((pos2.x-pos1.x)**2)+((pos2.y-pos1.y)**2)+((pos2.z-pos1.z)**2))
def accept(event, ax, fig, selector, all_selections):
    if event.key == "enter":
        print("Point selection OK!\nHere are the indices:")
        #all_selections.append(selector.ind.tolist()[0])
        all_selections.extend(selector.ind)
        print(all_selections)
        selector.disconnect()
        ax.set_title("(Selection OK, please close the window)")
        fig.canvas.draw()


title__ =""
calculatedistancetocounteratom =""

root = tk.Tk()
root.withdraw()
POSCARpath = filedialog.askopenfilename(title="Select the POSCAR file")
OUTCARpath = filedialog.askopenfilename(title="Select the OUTCAR file")
atoms = []
quantities = []
history = []
with open(POSCARpath, 'r') as POSCARfile:
    for line in POSCARfile:
        if line.strip() == "Selective dynamics":
            atoms = history[-2].split()
            quantities = [int(x) for x in history[-1].split()]
            break
        history.append(line.strip())
print("Atoms: " + str(atoms))
print("Quantities: " + str(quantities))
atomsymbolfulllist = []
for index_ in range(len(atoms)):
    for j in range(quantities[index_]):
        atomsymbolfulllist.append(atoms[index_])
print("Reading POSCAR file")

avector=None
bvector=None
cvector=None

with open(POSCARpath,"r") as inposfile:
    iscoordinate=False
    linecount=0
    for line in inposfile:
        linecount+=1
        if linecount==3:
            x, y, z = tuple(line.strip().split())
            x = float(x)
            y = float(y)
            z = float(z)
            avector=vpythvector(x, y, z)#np.array([float(k) for k in line.strip().split()])#/10
            print("a= ",avector,"")
            continue
        elif linecount==4:
            x, y, z = tuple(line.strip().split())
            x = float(x)
            y = float(y)
            z = float(z)
            bvector = vpythvector(x, y, z) #np.array([float(k) for k in line.strip().split()])#/10
            print("b= ", bvector,"")
            continue
        elif linecount==5:
            x, y, z = tuple(line.strip().split())
            x = float(x)
            y = float(y)
            z = float(z)
            cvector = vpythvector(x, y, z) #np.array([float(k) for k in line.strip().split()])#/10
            print("c= ", cvector,"")
            continue
        elif linecount>5:
            break

a=avector.mag
b=bvector.mag
c=cvector.mag
print("a length",a,"")
print("b length",b,"")
print("c length",c,"")
root.destroy()

wasposition = False
wasdotted = False
currentindex = -1
allatompos={}
molecule=[]
allatompos_firststep_OUTCARfile = None
time=0
firstatomsymbolnoid,secondatomsymbolnoid="",""

with open(OUTCARpath, 'r') as infile:
    for line in infile:
        if "POSITION" in line and "TOTAL-FORCE" in line:
            wasposition=True
            continue
        if "----------------------" in line and wasposition:
            wasposition=False
            wasdotted=True
            continue
        if wasdotted and not "------------------" in line:
            currentindex+=1
            pos_and_forces = line.strip().split()
            location_ = vpythvector(*[float(pos_and_forces[0]), float(pos_and_forces[1]), float(pos_and_forces[2])])# vpythvector(*(np.dot(np.array([float(pos_and_forces[0]), float(pos_and_forces[1]), float(pos_and_forces[2])]),latticevectorarray)))
            currentsymbol=atomsymbolfulllist[currentindex]+str(currentindex)
            allatompos[currentsymbol]=location_
        elif "-----------------------" in line and wasdotted:
            time+=1
            if time==1:
                allatompos_firststep_OUTCARfile=deepcopy(allatompos)
                for atomid,atomposlocation in allatompos.items():
                    molecule.append(atom(atomid, ''.join([i for i in atomid if not i.isdigit()]),atomposlocation))
                break
            currentindex=-1
            wasdotted=False
            continue
print(len(molecule))

# class SelectFromCollection(object):
#     def __init__(self, ax, collection, alpha_other=0.3):
#         self.canvas = ax.figure.canvas
#         self.collection = collection
#         self.alpha_other = alpha_other
#         self.xys = collection.get_offsets()
#         self.Npts = len(self.xys)
#         self.fc = collection.get_facecolors()
#         if len(self.fc) == 0:
#             raise ValueError('Collection must have a facecolor')
#         elif len(self.fc) == 1:
#             self.fc = np.tile(self.fc, (self.Npts, 1))
#         self.lasso = LassoSelector(ax, onselect=self.onselect)
#         self.ind = []
#         self.canvas.get_tk_widget().focus_set()
#     def onselect(self, verts):
#         path = Path(verts)
#         self.ind = np.nonzero(path.contains_points(self.xys))[0]
#         self.fc[:, -1] = self.alpha_other
#         self.fc[self.ind, -1] = 1
#         self.collection.set_facecolors(self.fc)
#         self.canvas.get_tk_widget().focus_set()
#         self.canvas.draw_idle()
#     def disconnect(self):
#         self.lasso.disconnect_events()
#         self.fc[:, -1] = 1
#         self.collection.set_facecolors(self.fc)
#         self.canvas.draw_idle()
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

# def accept(event):
#     global timesteps_selected
#     if event.key == "enter":
#         print("Points selected OK!")
#         timesteps_selected = selector.xys[
#             selector.ind]
#         print(timesteps_selected)
#         selector.disconnect()
#         plt.clf()
#         plt.close()
has_been_selected=False
ids_of_selectedatoms=[]
while not has_been_selected:
    ax, fig, pts = draw(molecule, set(atomsymbolfulllist), radiuslist, element_colors)
    fig.canvas.get_tk_widget().focus_set()
    selector = SelectFromCollection(ax, pts)
    fig.canvas.mpl_connect("key_press_event", lambda event: accept(event, ax, fig, selector, ids_of_selectedatoms))
    ax.set_title("Select SPE and press enter.\nOtherwise, close the window.")
    plt.show()
    plt.close()

    if confirm(text="Did you manage to select the atoms?", buttons=["Yes", "No"]) == "No":
        ids_of_selectedatoms = []
        continue
    else:
        ids_of_selectedatoms = selector.ind
        has_been_selected = True

    # orientation_option=confirm(text="Select the appropriate axis view.\nIf this is not the right choice, then press cancel in a later dialog to return here.",title="Select the axis orientation",buttons=["XY","YZ","XZ"])
    # data_firstaxis=[]
    # data_secondaxis=[]
    # sizes=[]
    # for atomi in molecule:
    #     if orientation_option == "XY":
    #         data_firstaxis.append(atomi.pos.x)
    #         data_secondaxis.append(atomi.pos.y)
    #         sizes.append(radiuslist[atomi.symbol]*10)
    #     elif orientation_option == "YZ":
    #         data_firstaxis.append(atomi.pos.z)
    #         data_secondaxis.append(atomi.pos.y)
    #         sizes.append(radiuslist[atomi.symbol] * 10)
    #     elif orientation_option == "XZ":
    #         data_firstaxis.append(atomi.pos.x)
    #         data_secondaxis.append(atomi.pos.z)
    #         sizes.append(radiuslist[atomi.symbol] * 10)
    #     else:
    #         raise NotImplementedError
    # ax = plt.gca()
    # fig=plt.gcf()
    # pts = ax.scatter(data_firstaxis, data_secondaxis, s=sizes)
    # fig.canvas.get_tk_widget().focus_set()
    # selector = SelectFromCollection(ax, pts)
    # def accept(event):
    #     global ids_of_selectedatoms
    #     if event.key == "enter":
    #         print("Point selection OK!\nHere are the indices:")
    #         ids_of_selectedatoms=selector.ind
    #         print(ids_of_selectedatoms)
    #         selector.disconnect()
    #         ax.set_title("(Selection OK, please close the window)")
    #         fig.canvas.draw()
    # fig.canvas.mpl_connect("key_press_event", accept)
    # ax.set_title("Select the molecule and press enter to accept selected points. Otherwise, close the window.")
    # plt.show()
    # plt.close()
    # if confirm(text="Did you manage to select the molecule?",title="Did you manage to select the molecule?",buttons=["Yes","No"]) == "Yes":
    #     has_been_selected=True
electrodeatomlist=[]
onlymoleculeatomlist=[]
count=-1
for atomindex in range(len(molecule)):
    count+=1
    temppos = molecule[atomindex].pos
    if count in ids_of_selectedatoms:
        onlymoleculeatomlist.append(molecule[atomindex])
    else:
        molecule[atomindex].iselectrode = True
        electrodeatomlist.append(molecule[atomindex])
electrodeatomids=[j.id for j in electrodeatomlist]
max_x=0
min_x=1e99
max_y=0
min_y=1e99
max_z=0
min_z=1e99
for electrodeatom in electrodeatomlist:
    if electrodeatom.pos.x>max_x:
        max_x=electrodeatom.pos.x
    if electrodeatom.pos.x<min_x:
        min_x=electrodeatom.pos.x
    if electrodeatom.pos.y>max_y:
        max_y=electrodeatom.pos.y
    if electrodeatom.pos.y<min_y:
        min_y=electrodeatom.pos.y
    if electrodeatom.pos.z>max_z:
        max_z=electrodeatom.pos.z
    if electrodeatom.pos.z<min_z:
        min_z=electrodeatom.pos.z
atom_closestconnections={}
for atom1 in onlymoleculeatomlist:
    atom_closestconnections[atom1.id]=[]
    for atom2 in onlymoleculeatomlist:
        if not atom1.id==atom2.id:
            atom1_symbol=atom1.symbol
            atom2_symbol=atom2.symbol
            atom1_radius=radiuslist[atom1_symbol]
            atom2_radius=radiuslist[atom2_symbol]
            #print(euclidiandistance_vectors(atom1.pos,atom2.pos))
            if euclidiandistance_vectors(atom1.pos,atom2.pos)<1.7:
                print(atom1.id,atom2.id)
                atom_closestconnections[atom1.id].append((atom2.id,atom2.symbol, euclidiandistance_vectors(atom1.pos,atom2.pos)))

atom_closestconnections = {k:v for k,v in atom_closestconnections.items() if v}
sorted_atomclosestconnections = {x:sorted(atom_closestconnections[x]) for x in atom_closestconnections.keys()}
overview=[]
final_atombondingpairslist=[]

for atom1id,bondinfolist in sorted_atomclosestconnections.items():
    if "H" in atom1id:
        if bondinfolist[0][0]+"-"+atom1id in overview or atom1id+"-"+bondinfolist[0][0] in overview:
            pass
        else:
            overview.append(bondinfolist[0][0]+"-"+atom1id)
            final_atombondingpairslist.append((atom1id,bondinfolist[0][0]))
    elif "C" in atom1id or "N" in atom1id or "O" in atom1id:
        if ("C" in atom1id and len(bondinfolist)>4) or ("O" in atom1id and len(bondinfolist)>2) or ("N" in atom1id and len(bondinfolist)>3):
            pass
        for bondpartner in bondinfolist:
            if bondpartner[0] + "-" + atom1id in overview or atom1id + "-" + bondpartner[0] in overview:
                pass
            else:
                overview.append(bondpartner[0] + "-" + atom1id)
                final_atombondingpairslist.append((atom1id,bondpartner[0]))
    else:
        #not implemented
        pass
print("There are ",len(final_atombondingpairslist)," bonds")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for atomi in onlymoleculeatomlist:
    ax.scatter3D(atomi.pos.x, atomi.pos.y, atomi.pos.z, color=colourdict[atomi.symbol],s=radiuslist[atomi.symbol]*100)
for bondatom1,bondatom2 in final_atombondingpairslist:
    firstatompos=[]
    secondatompos=[]
    for atom____ in onlymoleculeatomlist:
        if atom____.id==bondatom1:
            firstatompos=atom____.pos
        elif atom____.id==bondatom2:
            secondatompos=atom____.pos
    ax.plot([firstatompos.x,secondatompos.x],[firstatompos.y,secondatompos.y],[firstatompos.z,secondatompos.z],linewidth=6)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.title("View the bonded structure and check whether any misbonding occurred")

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize / 2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
axisEqual3D(ax)

plt.show()
plt.close()




###################################### PROCESS BOND LENGTH VARIATION IN OUTCAR FILE ###################################
print("\nReading through all simulation frames in OUTCAR... this may take up to five minutes...\n")
wasposition = False
wasdotted = False
firstinitialisationpassed = False
blenderobjectatomlist = {}
every_n_timesteps=1
currentindex = -1
allatompos={}
bondlengthvariationovertimedict={}
angletoelectrodeovertimedict={}
middleofbondyvaluesovertimedict={}
for bondpair________ in final_atombondingpairslist:
    bondlengthvariationovertimedict[bondpair________]=[]
    angletoelectrodeovertimedict[bondpair________]=[]
    middleofbondyvaluesovertimedict[bondpair________]=[]
proximitytoelectrodeatomsovertimedict={}
for moleculeatom in onlymoleculeatomlist:
    proximitytoelectrodeatomsovertimedict[moleculeatom.id]=[]
allatompos_firststep_OUTCARfile = None
negativeYaxisvector=vpythvector(0,-1,0)
time=0
firstatomsymbolnoid,secondatomsymbolnoid="",""
with open(OUTCARpath, 'r') as infile:
    for line in infile:
        if "POSITION" in line and "TOTAL-FORCE" in line:
            wasposition=True
            continue
        if "----------------------" in line and wasposition:
            wasposition=False
            wasdotted=True
            continue
        if wasdotted and not "------------------" in line:
            currentindex+=1
            pos_and_forces = line.strip().split()
            if len(pos_and_forces) >= 3 and currentindex in ids_of_selectedatoms:
                location_ = vpythvector(*[float(pos_and_forces[0]), float(pos_and_forces[1]), float(pos_and_forces[2])])#coords in outcar file are in cartesian (unlike poscar which is in direct)
                currentsymbol=atomsymbolfulllist[currentindex]+str(currentindex)
                allatompos[currentsymbol]=location_

        elif "-----------------------" in line and wasdotted:
            time+=1
            if time==1:
                allatompos_firststep_OUTCARfile=deepcopy(allatompos)
            for bondpairatom1,bondpairatom2 in final_atombondingpairslist:
                if time>=5000 and time<=15000:
                    atom1pos=allatompos[bondpairatom1]
                    atom2pos=allatompos[bondpairatom2]
                    # Orthogonal cell case
                    if abs(dot(avector, bvector)) <= 0.0001:
                        if abs(atom1pos.x - atom2pos.x) > a * 0.7:
                            if atom1pos.x>atom2pos.x:
                                atom1pos.x-=a
                            else:
                                atom2pos.x-=a
                        if abs(atom1pos.y - atom2pos.y) > b * 0.7:
                            if atom1pos.y>atom2pos.y:
                                atom1pos.y-=b
                            else:
                                atom2pos.y-=b
                        if abs(atom1pos.z - atom2pos.z) > c * 0.8:
                            if atom1pos.z>atom2pos.z:
                                atom1pos.z-=c
                            else:
                                atom2pos.z-=c

                    # Non-orthogonal cell cases
                    else:
                        abs_diff_x = abs(atom1pos.x - atom2pos.x)
                        abs_diff_y = abs(atom1pos.y - atom2pos.y)
                        if abs_diff_x > a * 0.7 and abs_diff_y < abs_diff_x:
                            if atom1pos.x > atom2pos.x:
                                atom1pos -= avector
                            else:
                                atom2pos -= avector
                        if abs_diff_y > b * 0.7 and abs_diff_x < abs_diff_y:
                            if atom1pos.y > atom2pos.y:
                                atom1pos -= bvector
                            else:
                                atom2pos -= bvector


                    atom1pos2pos_subtractionvector=atom2pos-atom1pos
                    bondlengthvariationovertimedict[(bondpairatom1,bondpairatom2)].append(atom1pos2pos_subtractionvector.mag)
            currentindex=-1
            wasdotted=False
            continue
print("Finished reading OUTCAR file")
print("Now at the RBC calculation")

###################################### PLOT STATISTICAL INFORMATION ####################################
print("\n££££££££££££££££££££££££££££££ Bond stretching result output data £££££££££££££££££££££££££££££££££\n££££££ See the output.csv file")
print('What is the cathode material?')
Cathode_Mat = input()
print('What is the name of the molecule?')
Name_Mol = input()
with open(f"output_{Cathode_Mat}_{Name_Mol}.csv","w") as openf:
    writer = csv.writer(openf)
    print("Bond pair type\tMode percentual bond length change distribution (ModeRBC) (%)\tMaximum of the percentual bond length change distribution (MaxRBC) (%)\tMinimum of the percentual bond length change distribution (MinRBC) (%)")
    writer.writerow("Bond pair type#Mode percentual bond length change distribution (ModeRBC) (%)#Maximum of the percentual bond length change distribution (MaxRBC) (%)#Minimum of the percentual bond length change distribution (MinRBC) (%)".split("#"))
    averages=[]
    stdevs=[]
    maxes=[]
    mins=[]
    for bondpairtuple, instantaneousbondlengthvariation in bondlengthvariationovertimedict.items():
        firstatompos = []
        secondatompos = []
        for atom_ in onlymoleculeatomlist:
            if atom_.id == bondpairtuple[0]:
                firstatompos = atom_.pos
                continue
            elif atom_.id == bondpairtuple[1]:
                secondatompos = atom_.pos
                continue
        bondlengthinvacuum = euclidiandistance_vectors(firstatompos, secondatompos)
        normalizedbondlengthvariation= [((l/bondlengthinvacuum)-1)*100 for l in instantaneousbondlengthvariation]
        nbins = 100
        n, bins, patches = plt.hist(normalizedbondlengthvariation, density=True,bins=nbins)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        APBLI = round(bin_centers[n.argmax()], 2)
        averagebondlengthincrease=APBLI
        averages.append(averagebondlengthincrease)
        maxofbondlengthincrease = round(np.nanmax(normalizedbondlengthvariation), 2)
        maxes.append(maxofbondlengthincrease)
        minofbondlengthincrease = round(np.nanmin(normalizedbondlengthvariation), 2)
        mins.append(minofbondlengthincrease)
        print(bondpairtuple[0]+"-"+bondpairtuple[1]+"\t"+str(averagebondlengthincrease)+"\t"+str(maxofbondlengthincrease)+"\t"+str(minofbondlengthincrease))
        writer.writerow([bondpairtuple[0]+"-"+bondpairtuple[1],str(averagebondlengthincrease),str(maxofbondlengthincrease),str(minofbondlengthincrease)])
print("\n££££££££££££££££££££££££££££££ end £££££££££££££££££££££££££££££££££")

###################################### SEPARATE BOND PLOTTING out of interest ###################################
#Plot each bond property separately
count=-1
for bondatompair in final_atombondingpairslist:
    print(bondatompair)
    firstatomid,secondatomid=bondatompair
    count+=1
    firstatompos = []
    secondatompos = []
    for atom_ in onlymoleculeatomlist:
        if atom_.id == firstatomid:
            firstatompos = atom_.pos
            continue
        elif atom_.id == secondatomid:
            secondatompos = atom_.pos
            continue
    bondlengthinvacuum = euclidiandistance_vectors(firstatompos, secondatompos)
    increaseinbondlengthrefinitialvacuumbondlength = [((l / bondlengthinvacuum) - 1) * 100 for l in
                                                      bondlengthvariationovertimedict[bondatompair]]

    # Create a figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # First subplot
    ax1.plot(increaseinbondlengthrefinitialvacuumbondlength,
             label="Percentual increase in bond length relative to vacuum", color=(226 / 255, 58 / 255, 75 / 255))
    cm = plt.cm.get_cmap("hsv")
    nbins = 100
    n, bins, patches = ax2.hist(increaseinbondlengthrefinitialvacuumbondlength, density=True, bins=nbins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    counter = -1
    for i, p in enumerate(patches):
        counter += 1
        plt.setp(p, 'facecolor', cm((bin_centers[counter] + 7) / (7 * 2)))
    ax2.set_ylabel('Probability', fontsize=10, fontweight='bold')
    ax2.set_xlabel('Percentual increase in bond length ref. vacuum bond length', fontsize=10, fontweight='bold')
    ax2.set_title("Distribution for " + str(bondatompair[0]) + " and " + str(bondatompair[1]))
    ax2.set_xlim(min(increaseinbondlengthrefinitialvacuumbondlength),
                 max(increaseinbondlengthrefinitialvacuumbondlength))

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    plt.savefig(Cathode_Mat + "_" + Name_Mol + "_" + str(bondatompair[0]) + " " + str(bondatompair[1]) + '.png',bbox_inches='tight', dpi=500)

    #plt.show()
    plt.close()
print("Finished :D")