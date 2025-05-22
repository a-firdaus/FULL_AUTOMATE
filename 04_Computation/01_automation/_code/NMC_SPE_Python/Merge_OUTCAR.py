import sys
import time
import os
import re


# This v2 version skip the first lines of the OUTCAR files that are not the first ones, i.e. second one, third one, etc.
# The Iteration number are modified in the new OUTCAR so that there is no interruption when going from one file to the next one

def concatenate_outcars(file_list, output_file):
    start_time = time.time()

    def print_loading_bar(iteration, total, bar_length=50):
        progress = iteration / total
        block = int(bar_length * progress)
        bar = "#" * block + "-" * (bar_length - block)
        sys.stdout.write(f"\rProgress: [{bar}] {progress:.2%}")
        sys.stdout.flush()

    # Calculate the total number of lines
    total_lines = 15000
    # for file in file_list:
    #     with open(file, 'r') as f:
    #         total_lines += sum(1 for _ in f)

    lines_processed = 0
    firstLine = 'First call to EWALD'
    canStartWritting = False

    with open(output_file, 'w') as output:
        iterList = []
        counter = 1
        for i, file in enumerate(file_list):
            with open(file, 'r') as f:
                print(f'\nFile number {i + 1}', '\nFile name :', file)
                #print('\n',file)
                lines = f.readlines()
                if i == 0:
                    for j, line in enumerate(lines):
                        output.write(line)
                        # lines_processed += 1
                        # print_loading_bar(lines_processed, total_lines)
                        if "Iteration" in line:
                            if len(iterList) == 0:
                                iterList.append(int(line.split()[2][:-1]))
                                lines_processed += 1
                                print_loading_bar(lines_processed, total_lines)
                            else:
                                newIter = int(line.split()[2][:-1])
                                if newIter < int(iterList[-1]):
                                    iterList.append(newIter)
                                    lines_processed += 1
                                    print_loading_bar(lines_processed, total_lines)
                                elif newIter > int(iterList[-1]):
                                    iterList.pop(-1)
                                    iterList.append(newIter)
                                    lines_processed += 1
                                    print_loading_bar(lines_processed, total_lines)
                else:
                    lineNumber = 0
                    for j, line in enumerate(lines):
                        start_iteration = sum(iterList)
                        if firstLine in line:
                            lineNumber = j
                            #print(f'<<First call to EWALD>> was found in the current line (#{j})')
                        elif lineNumber != 0 and j > lineNumber+2:
                            #if j == lineNumber+3:
                                #print(f'Two lines after the triggering line was skipped. We are now on line {j}')
                            if "Iteration" in line:
                                newIter = int(line.split()[2][:-1])
                                if newIter < int(iterList[-1]):
                                    iterList.append(newIter)
                                    output.write(re.sub(r'\d+', str(sum(iterList)), line, count=1))
                                    lines_processed += 1
                                    print_loading_bar(lines_processed, total_lines)
                                elif newIter > int(iterList[-1]):
                                    iterList.pop(-1)
                                    iterList.append(newIter)
                                    output.write(re.sub(r'\d+', str(sum(iterList)), line, count=1))
                                    lines_processed += 1
                                    print_loading_bar(lines_processed, total_lines)
                            else:
                                output.write(line)

            output.write('\n')
                #output.writelines(lines)
        print("\nNumber of iterations (=steps) is :",sum(iterList))
    end_time = time.time()
    print(f"Total time taken: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    # Provide the paths to the OUTCAR files from your runs
    outcar_file_list = [os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR_1'),
                        os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR_2'),
                        os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR_3'),
                        os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR_4'),
                        os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR_5'),
                        os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR_6'),
                        ]


    # Provide the path where you want to save the concatenated OUTCAR-like file
    output_contcar_file = os.path.join(os.path.expanduser('~'), 'Desktop', 'SPE_NMC', 'With_EF_simulations', 'LiNMC811', 'Alcohol', 'OUTCAR')
    concatenate_outcars(outcar_file_list, output_contcar_file)