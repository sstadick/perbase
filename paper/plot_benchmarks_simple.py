#!/usr/bin/env python3
"""
Generate benchmark plots for the JOSS paper - simplified version
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load data
data_dir = Path("data/outputs")
base_depth_df = pd.read_csv(data_dir / "bench_results_base_depth.csv")
mate_fix_df = pd.read_csv(data_dir / "bench_results_base_depth_mate_fix.csv")

# Extract tool names
base_depth_df['tool'] = base_depth_df['command'].apply(lambda x: x.split()[0])
mate_fix_df['tool'] = mate_fix_df['command'].apply(lambda x: x.split()[0])

# Create output directory
output_dir = Path("outputs")
output_dir.mkdir(exist_ok=True)

# Create a simple comparison plot
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Prepare data
tools = ['perbase', 'sambamba']
modes = ['Standard', 'Mate-fix']

perbase_std = base_depth_df[base_depth_df['tool'] == 'perbase']['mean'].values[0]
sambamba_std = base_depth_df[base_depth_df['tool'] == 'sambamba']['mean'].values[0]
perbase_mf = mate_fix_df[mate_fix_df['tool'] == 'perbase']['mean'].values[0]
sambamba_mf = mate_fix_df[mate_fix_df['tool'] == 'sambamba']['mean'].values[0]

x = np.arange(len(tools))
width = 0.35

bars1 = ax.bar(x - width/2, [perbase_std, sambamba_std], width, label='Standard')
bars2 = ax.bar(x + width/2, [perbase_mf, sambamba_mf], width, label='Mate-fix')

ax.set_ylabel('Runtime (ms)')
ax.set_xlabel('Tool')
ax.set_title('perbase vs sambamba: Runtime Performance')
ax.set_xticks(x)
ax.set_xticklabels(tools)
ax.legend()

# Add speedup text
speedup_std = sambamba_std / perbase_std
speedup_mf = sambamba_mf / perbase_mf

ax.text(0.5, 0.95, f'Standard mode: perbase is {speedup_std:.1f}x faster', 
        ha='center', va='top', transform=ax.transAxes,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Save plot
plt.savefig(output_dir / 'benchmark_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'benchmark_comparison.pdf', bbox_inches='tight')

print("Plot generated successfully!")
print(f"Standard mode: perbase is {speedup_std:.1f}x faster than sambamba")
print(f"Mate-fix mode: perbase is {speedup_mf:.1f}x faster than sambamba")