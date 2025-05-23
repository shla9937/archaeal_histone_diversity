#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load your CSV
file = "rmsds.csv"  # Replace with your CSV file path
df = pd.read_csv(file, index_col=0)

# Convert strings to numeric, ignore errors (i.e. empty cells)
df = df.apply(pd.to_numeric, errors='coerce')

# Calculate mean and SEM (Standard Error of the Mean)
means = df.mean(axis=1)
sems = df.sem(axis=1)

# Custom names and colors
custom_names = {
    "control_euk": "Eukaryotic",
    "control_arc": "Archaeal",
    "single_1": "Single 1",
    "single_2": "Single 2",
    "single_3": "Single 3",
    "multi_1": "Multi 1",
    "multi_1_1": "Multi 1: A",
    "multi_1_2": "Multi 1: B",
    "combo_12": "Combo 1&2",
    "combo_12_1": "Combo 1&2: 1",
    "combo_12_2": "Combo 1&2: 2",
    "combo_34_handmade": "Combo 3&4",
    "combo_34_3": "Combo 3&4: 3"
}

custom_colors = {
    "control_euk": "000000",
    "control_arc": "#000000",
    "single_1": "#56B4E9",
    "single_2": "#D55E00",
    "single_3": "#E69F00",
    "multi_1": "#0072B2",
    "multi_1_1": "#4095C5",
    "multi_1_2": "#0072B2",
    "combo_12": "#AAAAAA",
    "combo_12_1": "#56B4E9",
    "combo_12_2": "#D55E00",
    "combo_34_handmade": "#AAAAAA",
    "combo_34_3": "#E69F00"
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
ax.set_title("RMSD from first frame to average last 10ns", fontsize=6)
plt.yticks(fontsize=6)
plt.xticks(fontsize=6)
plt.tight_layout()
plt.savefig('combined_finals/outputs/nuc_rmsds.png', dpi=300)
plt.show()
