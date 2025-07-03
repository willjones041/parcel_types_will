import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  
import os
import shutil

# Prompt for number of datasets
num_datasets = int(input("How many datasets do you want to plot? (2 or more): "))

# Store all data in lists
all_times = []
all_zpos = []
all_qr = []
all_nr = []
all_vterm = []
all_names = []

for i in range(num_datasets):
    filename = input(f"Enter the path to CSV file #{i+1}: ")
    print(f"Loading data from {filename}")
    df = pd.read_csv(filename)
    # Assume time is col 0, zpos is col 1, qr col 8, nr col 9, vterm col 10
    time = df.iloc[:, 0]
    zpos = df.iloc[:, 1]
    qr = df.iloc[:, 8]
    nr = df.iloc[:, 9]
    vterm = df.iloc[:, 10]
    all_times.append(time)
    all_zpos.append(zpos)
    all_qr.append(qr)
    all_nr.append(nr)
    all_vterm.append(vterm)
    name = input(f"Enter a label for dataset #{i+1}: ")
    all_names.append(name)

foldername = input("Enter a folder name to save the plots: ")
if os.path.exists(foldername):
    shutil.rmtree(foldername)
os.makedirs(foldername, exist_ok=True)

sns.set()

# Height vs time
plt.figure(figsize=(10, 6))
colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k']
for i in range(num_datasets):
    plt.plot(all_times[i], all_zpos[i], marker='o', color=colors[i % len(colors)], markersize=1)
plt.legend(all_names, fontsize=12)
plt.xlabel('Time (s)', fontsize=24)
plt.ylabel('Height (m)', fontsize=24)
plt.tick_params(axis='both', labelsize=20)
plt.title("Height vs time", fontsize=24)
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(foldername, f"{foldername}_Height_vs_time.svg"))
plt.show()

# qr vs height
plt.figure(figsize=(10, 6))
for i in range(num_datasets):
    plt.plot(all_zpos[i], all_qr[i], marker='o', color=colors[i % len(colors)], markersize=1)
plt.legend(all_names, fontsize=12)
plt.xlabel('Height (m)', fontsize=24)
plt.ylabel('Rainwater mixing ratio (kg kg$^{-1}$)', fontsize=24)
plt.tick_params(axis='both', labelsize=20)
plt.title("Rainwater mixing ratio vs height", fontsize=24)
plt.gca().invert_xaxis()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(foldername, f"{foldername}_qr_vs_height.svg"))
plt.show()

# nr vs height
plt.figure(figsize=(10, 6))
for i in range(num_datasets):
    plt.plot(all_zpos[i], all_nr[i], marker='o', color=colors[i % len(colors)], markersize=1)
plt.legend(all_names, fontsize=12)
plt.xlabel('Height (m)', fontsize=24)
plt.ylabel('Rainwater number concentration (m$^{{-3}}$)', fontsize=16)
plt.tick_params(axis='both', labelsize=20)
plt.title("Rainwater number concentration vs height", fontsize=24)
plt.gca().invert_xaxis()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(foldername, f"{foldername}_nr_vs_height.svg"))
plt.show()

# vterm vs height
plt.figure(figsize=(10, 6))
for i in range(num_datasets):
    plt.plot(all_zpos[i], all_vterm[i], marker='o', color=colors[i % len(colors)], markersize=1)
plt.legend(all_names, fontsize=12)
plt.xlabel('Height (m)', fontsize=24)
plt.ylabel('Terminal velocity (ms$^{-1}$)', fontsize=24)
plt.tick_params(axis='both', labelsize=20)
plt.title("Terminal velocity vs height", fontsize=24)
plt.gca().invert_xaxis()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(foldername, f"{foldername}_Terminal_velocity_vs_height.svg"))
plt.show()