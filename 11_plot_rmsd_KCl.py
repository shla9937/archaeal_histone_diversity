#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load your CSV
file = "rmsds_KCl.csv"  # Replace with your CSV file path
df = pd.read_csv(file, index_col=0)

# Convert strings to numeric, ignore errors (i.e. empty cells)
df = df.apply(pd.to_numeric, errors='coerce')

# Calculate mean and SEM (Standard Error of the Mean)
means = df.mean(axis=1)
sems = df.sem(axis=1)

# Custom names and colors
custom_names = {
    "single_3": "Single 3",
    "single_3_2M_KCl": "Single 3 - 2M KCl",
    "multi_1_1": "Multi 1: A",
    "multi_1_1_2M_KCl": "Multi 1: A - 2M KCl",
    "combo_34_handmade": "Combo 3&4",
    "combo_34_handmade_2M_KCl": "Combo 3&4 - 2M KCl",
    "combo_34_3": "Combo 3&4: 3",
    "combo_34_3_2M_KCl": "Combo 3&4: 3 - 2M KCl"
}

custom_colors = {
    "single_3": "#444444",
    "single_3_2M_KCl": "#E69F00",
    "multi_1_1": "#444444",
    "multi_1_1_2M_KCl": "#4095C5",
    "combo_34_handmade": "#444444",
    "combo_34_handmade_2M_KCl": "#AAAAAA",
    "combo_34_3": "#444444",
    "combo_34_3_2M_KCl": "#E69F00"
}

# Plotting
fig, ax = plt.subplots(figsize=(3, 2))    
plt.rcParams.update({'font.size': 6})
ax.spines[['right', 'top']].set_visible(False)

x = np.arange(len(means))
bar_colors = [custom_colors.get(i, "#000000") for i in means.index]
bar_labels = [custom_names.get(i, i) for i in means.index]

bars = ax.bar(x, means, yerr=sems, capsize=5, color=bar_colors)

ax.set_xticks(x)
ax.set_xticklabels(bar_labels, rotation=90, ha='right')
ax.set_ylabel("RMSD (Å)", fontsize=6)
ax.set_title("Average RMSD of DNA from first frame to last 10ns", fontsize=6)
plt.yticks(fontsize=6)
plt.xticks(fontsize=6)
plt.tight_layout()
plt.savefig('combined_finals/outputs/nuc_rmsds_KCl.png', dpi=300)
plt.show()
