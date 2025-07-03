import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data from file
with open('atmos_profile.dat', 'r') as f:
    header = f.readline().strip().replace('#','').split()

data = np.loadtxt('atmos_profile.dat', skiprows=1)

# First column is x-axis (e.g., height or pressure)
x = data[:, 0]

# Number of columns (excluding the first)
num_cols = data.shape[1]

# Create subplots for each column (except the first)
rows = 3
cols = 2
fig, axes = plt.subplots(rows, cols, figsize=(12, 12), sharex=True)
axes = axes.flatten()

sns.set(style="whitegrid", context="notebook", palette="deep")

for i in range(1, min(num_cols, rows*cols)+1):
    sns.lineplot(x=x, y=data[:, i], ax=axes[i - 1])
    axes[i - 1].set_ylabel(header[i])
    axes[i - 1].grid(True)

# Hide any unused subplots
for j in range(num_cols, rows*cols+1):
    axes[j-1].set_visible(False)

for ax in axes:
    ax.set_xlabel(header[0])

plt.tight_layout()
plt.show()