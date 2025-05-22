import warnings
warnings.filterwarnings("ignore")
from math import sqrt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from vpython import vector as vpythvector
import tkinter as tk
from tkinter import filedialog
from molatomclasses import atom
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import os
from collections import Counter
from pyautogui import confirm
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

def euclidiandistance_vectors(pos1,pos2):
    return sqrt(((pos2.x-pos1.x)**2)+((pos2.y-pos1.y)**2)+((pos2.z-pos1.z)**2))

root = tk.Tk()
root.withdraw()
POSCARpath = filedialog.askopenfilename(title="Select the POSCAR file")
CONTCARpath = filedialog.askopenfilename(title="Select the CONTCAR file")

rewrite = input("Do you want to 're-write' the POSCAR (1) or CONTCAR (2)?\n")
while rewrite not in ["1", "2"]:
    print("Please enter a valid input! 1 for POSCAR or 2 for CONTCAR")
    rewrite = input()

file = CONTCARpath if rewrite=="2" else POSCARpath

directory_path = os.path.dirname(POSCARpath)
folder_name = os.path.basename(directory_path)
parent_directory = os.path.basename(os.path.dirname(directory_path))
location = parent_directory + "_" + folder_name[:3]
print("Location:", location)

atoms = []
quantities = []
history = []
nbr_line = -9
check = False
with open(POSCARpath, 'r') as POSCARfile:
    for line in POSCARfile:
        nbr_line += 1
        if line.strip() == "Selective dynamics":
            check = True
            atoms = history[-2].split()
            quantities = [int(x) for x in history[-1].split()]
        if not check:
            history.append(line.strip())
print("Number of lines in the POSCAR file : ", nbr_line)
print("Atoms: " + str(atoms))
print("Quantities: " + str(quantities))
atomsymbolfulllist = []
for index_ in range(len(atoms)):
    for j in range(quantities[index_]):
        atomsymbolfulllist.append(atoms[index_])

avector=None
bvector=None
cvector=None
a = None
b = None
c = None
latticevectorarray=None

currentindex = -1
allatompos={}
molecule=[]

with open(POSCARpath,"r") as inposfile:
    iscoordinate=False
    linecount=0
    for line in inposfile:
        linecount+=1
        if linecount==3:
            avector=np.array([float(k) for k in line.strip().split()])
            a = np.linalg.norm(avector)
            print("a= ",avector,"")
            continue
        elif linecount==4:
            bvector = np.array([float(k) for k in line.strip().split()])
            b = np.linalg.norm(bvector)
            print("b= ", bvector,"")
            continue
        elif linecount==5:
            cvector = np.array([float(k) for k in line.strip().split()])
            print("c= ", cvector,"")
            c = np.linalg.norm(cvector)
            latticevectorarray=np.stack([avector,bvector,cvector], axis=0)
            print("latticevectorarray=",latticevectorarray,"")
            continue
        elif linecount>9 and linecount<sum(quantities, 10):
            currentindex += 1
            pos = line.strip().split()
            location_ = vpythvector(*[float(pos[0])*a, float(pos[1])*b, float(pos[2])*c])
            currentsymbol = atomsymbolfulllist[currentindex] + str(currentindex)
            allatompos[currentsymbol] = location_
    for atomid, atomposlocation in allatompos.items():
        #print(atomid, ''.join([i for i in atomid if not i.isdigit()]))
        molecule.append(atom(atomid, ''.join([i for i in atomid if not i.isdigit()]), atomposlocation))
    #print(molecule)

root.destroy()

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

has_been_selected=False
ids_of_selectedatoms=[]
while not has_been_selected:
    orientation_option=confirm(text="Select the appropriate axis view.\nIf this is not the right choice, then press cancel in a later dialog to return here.",
                               title="Select the axis orientation",buttons=["XY","YZ","XZ"])
    data_firstaxis=[]
    data_secondaxis=[]
    sizes=[]
    for atomi in molecule:
        if orientation_option == "XY":
            data_firstaxis.append(atomi.pos.x)
            data_secondaxis.append(atomi.pos.y)
            sizes.append(radiuslist[atomi.symbol]*10)
        elif orientation_option == "YZ":
            data_firstaxis.append(atomi.pos.z)
            data_secondaxis.append(atomi.pos.y)
            sizes.append(radiuslist[atomi.symbol] * 10)
        elif orientation_option == "XZ":
            data_firstaxis.append(atomi.pos.x)
            data_secondaxis.append(atomi.pos.z)
            sizes.append(radiuslist[atomi.symbol] * 10)
        else:
            raise NotImplementedError
    ax = plt.gca()
    fig=plt.gcf()
    pts = ax.scatter(data_firstaxis, data_secondaxis, s=sizes)
    fig.canvas.get_tk_widget().focus_set()
    selector = SelectFromCollection(ax, pts)
    def accept(event):
        global ids_of_selectedatoms
        if event.key == "enter":
            print("Point selection OK!\nHere are the indices:")
            ids_of_selectedatoms=selector.ind
            print(ids_of_selectedatoms)
            selector.disconnect()
            ax.set_title("(Selection OK, please close the window)")
            fig.canvas.draw()
    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Select the molecule and press enter to accept selected points. Otherwise, close the window.")
    plt.show()
    plt.close()
    if confirm(text="Did you manage to select the molecule?",title="Did you manage to select the molecule?",buttons=["Yes","No"]) == "Yes":
        has_been_selected=True
electrodeatomlist=[]
onlymoleculeatomlist=[]
toBeDeletedList_symbol = []
count=-1
for atomindex in range(len(molecule)):
    count+=1
    temppos = molecule[atomindex].pos
    if count in ids_of_selectedatoms:
        onlymoleculeatomlist.append(molecule[atomindex])
        toBeDeletedList_symbol.append(molecule[atomindex].symbol)
    else:
        molecule[atomindex].iselectrode = True
        electrodeatomlist.append(molecule[atomindex])
print("Atoms selected to be deleted:", toBeDeletedList_symbol)
toBeDeleted_count = dict(Counter(toBeDeletedList_symbol))
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
            if euclidiandistance_vectors(atom1.pos,atom2.pos)<1.7:
                #print(atom1.id,atom2.id, ":", euclidiandistance_vectors(atom1.pos,atom2.pos))
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
#print("There are ",len(final_atombondingpairslist)," bonds")

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

if len(onlymoleculeatomlist) > len(molecule)/2:
    name = "_Molecule"
else:
    name = "_Electrode"
new_CONTCAR_path = location + name

with open(file, 'r') as original_file, open(new_CONTCAR_path, 'w') as new_file:
    counter  = 0
    total = 0
    for i, line in enumerate(original_file):
        if i == 0:
            new_file.write(new_CONTCAR_path + '\n')
        elif i > 0 and i < 5:
            new_file.write(line)  # Keep the first 5 lines as they are
        elif i == 5:
            elements_line = line.strip().split()
            new_elements_line = []
            for i, elem in enumerate(elements_line):
                if elem not in toBeDeleted_count.keys():
                    new_elements_line.append(elem)
                else:
                    q_init = int(quantities[i])
                    q_toDelete = int(toBeDeleted_count[elem])
                    if q_toDelete < q_init:
                        new_elements_line.append(elem)

            new_file.write(" ".join(new_elements_line) + "\n")  # Modify the 6th line
        elif i == 6:
            quantities_line = line.strip().split()
            new_quantities_line = []
            for i, q in enumerate(quantities_line):
                elem = atoms[i]
                if elem not in toBeDeleted_count.keys():
                    new_quantities_line.append(q)
                else:
                    temp = int(q)
                    if temp > int(toBeDeleted_count[elem]):
                        temp -= int(toBeDeleted_count[elem])
                        new_quantities_line.append(f"{temp}")

            new_file.write(" ".join(new_quantities_line) + "\n")  # Modify the 7th line
            total = sum([int(x) for x in new_quantities_line])
        elif i==7 or i==8:
            new_file.write(line)  # Keep lines 8 and 9 as they are
        elif i >= 9 and i<nbr_line + 9:
            if i - 9 not in [ids for ids in ids_of_selectedatoms]:
                new_file.write(line)
                counter += 1

print("Number of elements matches the number of lines:", total == counter)
print(f"New CONTCAR file created: {new_CONTCAR_path}")