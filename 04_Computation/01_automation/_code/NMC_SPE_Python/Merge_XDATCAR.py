import os
import tkinter as tk
from tkinter import filedialog, simpledialog

def read_xdatcar(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines
def gather_xdatcar_data(num_files):
    data = []

    # Use file dialog to select XDATCAR files
    root = tk.Tk()
    root.withdraw()

    remaining_files = num_files
    file_paths = []

    # Use the directory of the first selected file as the output directory
    output_directory = None

    while remaining_files > 0:
        file_path = filedialog.askopenfilename(
            initialdir="/path/to/xdatcar/files",
            title=f"Select XDATCAR file ({remaining_files} remaining)")
            #filetypes=[("XDATCAR files", "XDATCAR")]


        if file_path:
            file_paths.append(file_path)
            remaining_files -= 1

            # Extract the directory of the first selected file
            if output_directory is None:
                output_directory = os.path.dirname(file_path)

    # Specify the output file path in the first file's directory
    output_file = os.path.join(output_directory, "Combined_XDATCAR")

    # Loop through selected XDATCAR files
    for file_path in file_paths:
        filename = os.path.basename(file_path)
        xdatcar_data = read_xdatcar(file_path)
        data.append({
            'filename': filename,
            'data': xdatcar_data
        })

    # Write the gathered data to the specified output file
    with open(output_file, 'w') as output:
        modified_data = []
        last_iteration = 0
        counter = 0
        skip_counter = 0
        first_time = 0
        for entry in data:
            for line in entry['data']:
                if 'NMC' in line and first_time == 0:
                    first_time += 1
                    modified_data.append(line)
                else:
                    if 'NMC' in line:
                        skip_counter = 6
                        continue
                    if skip_counter > 0:
                        skip_counter -= 1
                        continue
                    elif "Direct configuration=" in line:
                        counter += 1
                        if last_iteration < int(line.split("=")[1]):
                            modified_data.append(line)
                            last_iteration = int(line.split("=")[1])
                        else:
                            last_iteration += 4
                            modified_line = f"Direct configuration= {last_iteration}\n"
                            modified_data.append(modified_line)
                    else:
                        modified_data.append(line)

        # Write the content of the XDATCAR file or any specific information you want
        output.writelines(modified_data)
        print(counter*4)


if __name__ == "__main__":
    output_file = "Combined_XDATCAR"

    # Use simpledialog to ask for the number of files
    root = tk.Tk()
    root.withdraw()
    num_files = simpledialog.askinteger("Number of Files", "How many files to gather?")

    gather_xdatcar_data(num_files)