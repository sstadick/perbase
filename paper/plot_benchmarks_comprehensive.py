#!/usr/bin/env python3
"""
Generate comprehensive benchmark plots for the JOSS paper including all mate-fix methods
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load data
data_dir = Path("data/outputs")
base_depth_df = pd.read_csv(data_dir / "bench_results_base_depth.csv")
mate_fix_df = pd.read_csv(data_dir / "bench_results_base_depth_mate_fix.csv")
mate_fix_new_df = pd.read_csv(data_dir / "bench_results_base_depth_mate_fix_new_v_old.csv")

# Create output directory
output_dir = Path("outputs")
output_dir.mkdir(exist_ok=True)

# Extract tool names and convert to minutes
def process_data(df):
    df = df.copy()
    df['tool'] = df['command'].apply(lambda x: x.split()[0])
    df['mean_minutes'] = df['mean'] / 60  # Convert to minutes
    return df

base_depth_df = process_data(base_depth_df)
mate_fix_df = process_data(mate_fix_df)
mate_fix_new_df = process_data(mate_fix_new_df)

# Get standard runtime data
perbase_std = base_depth_df[base_depth_df['tool'] == 'perbase']['mean_minutes'].values[0]
sambamba_std = base_depth_df[base_depth_df['tool'] == 'sambamba']['mean_minutes'].values[0]

# Get sambamba mate-fix runtime
sambamba_mf = mate_fix_new_df[mate_fix_new_df['command'].str.contains('sambamba')]['mean_minutes'].values[0]

# Extract mate-fix method data (excluding OLD version)
mate_fix_methods = []
mate_fix_runtimes = []

for _, row in mate_fix_new_df.iterrows():
    if 'perbase' in row['command'] and '30X OLD' not in row['command']:
        # Extract method name from command
        parts = row['command'].split()
        if len(parts) >= 4:
            method = ' '.join(parts[3:])  # Everything after "perbase base-depth mate-fix"
            method = method.replace('30X ', '')  # Remove "30X" prefix
            mate_fix_methods.append(method)
            mate_fix_runtimes.append(row['mean_minutes'])

# Create comprehensive comparison plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Simple comparison (Standard vs best mate-fix method)
best_mate_fix_runtime = min(mate_fix_runtimes)
best_mate_fix_method = mate_fix_methods[np.argmin(mate_fix_runtimes)]

modes = ['Standard', f'Mate-fix\n({best_mate_fix_method})']
x = np.arange(len(modes))
width = 0.35

bars1 = ax1.bar(x - width/2, [perbase_std, best_mate_fix_runtime], width, 
                label='perbase', color='#1f77b4')
bars2 = ax1.bar(x + width/2, [sambamba_std, sambamba_mf], width, 
                label='sambamba', color='#ff7f0e')

ax1.set_ylabel('Runtime (minutes)')
ax1.set_xlabel('Mode')
ax1.set_title('perbase vs sambamba: base-depth Performance')
ax1.set_xticks(x)
ax1.set_xticklabels(modes)
ax1.legend()

# Plot 2: All mate-fix methods comparison
# Sort methods by performance
sorted_indices = np.argsort(mate_fix_runtimes)
sorted_methods = [mate_fix_methods[i] for i in sorted_indices]
sorted_runtimes = [mate_fix_runtimes[i] for i in sorted_indices]

# Truncate method names for better display
display_methods = []
for method in sorted_methods:
    if len(method) > 20:
        # Abbreviate long method names
        if 'BaseQualMapQual' in method:
            method = method.replace('BaseQualMapQual', 'BQ→MQ→')
        elif 'MapQualBaseQual' in method:
            method = method.replace('MapQualBaseQual', 'MQ→BQ→')
        display_methods.append(method)
    else:
        display_methods.append(method)

y_pos = np.arange(len(sorted_methods))

# Horizontal bar chart for mate-fix methods
bars = ax2.barh(y_pos, sorted_runtimes, color='#2ca02c', alpha=0.7)
ax2.set_yticks(y_pos)
ax2.set_yticklabels(display_methods, fontsize=9)
ax2.set_xlabel('Runtime (minutes)')
ax2.set_title('perbase Mate-fix Method Performance Comparison')
ax2.grid(axis='x', alpha=0.3)

# Add runtime values on bars
for i, (bar, runtime) in enumerate(zip(bars, sorted_runtimes)):
    ax2.text(runtime + 0.5, bar.get_y() + bar.get_height()/2, 
             f'{runtime:.1f}m', ha='left', va='center', fontsize=8)

# Add reference line for standard mode
ax2.axvline(x=perbase_std, color='#1f77b4', linestyle='--', alpha=0.7, 
            label=f'Standard mode ({perbase_std:.1f}m)')
ax2.legend()

plt.tight_layout()
plt.savefig(output_dir / 'benchmark_comprehensive.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'benchmark_comprehensive.pdf', bbox_inches='tight')

# Create the simple version for the paper (same as before but with best mate-fix method)
fig2, ax = plt.subplots(1, 1, figsize=(8, 6))

x = np.arange(2)  # Standard and Mate-fix
width = 0.35

bars1 = ax.bar(x - width/2, [perbase_std, best_mate_fix_runtime], width, 
               label='perbase', color='#1f77b4')
bars2 = ax.bar(x + width/2, [sambamba_std, sambamba_mf], width, 
               label='sambamba', color='#ff7f0e')

ax.set_ylabel('Runtime (minutes)')
ax.set_xlabel('Mode')
ax.set_title('perbase vs sambamba: base-depth Performance')
ax.set_xticks(x)
ax.set_xticklabels(['Standard', 'Mate-fix'])
ax.legend()

plt.tight_layout()
plt.savefig(output_dir / 'benchmark_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'benchmark_comparison.pdf', bbox_inches='tight')
plt.close()

# Calculate speedups
speedup_std = sambamba_std / perbase_std
speedup_mf = sambamba_mf / best_mate_fix_runtime

print("Plots generated successfully!")
print("\nPerformance Summary:")
print(f"Standard mode:")
print(f"  perbase: {perbase_std:.1f} minutes")
print(f"  sambamba: {sambamba_std:.1f} minutes")
print(f"  Speedup: {speedup_std:.1f}x faster")
print(f"\nMate-fix mode (best method: {best_mate_fix_method}):")
print(f"  perbase: {best_mate_fix_runtime:.1f} minutes")
print(f"  sambamba: {sambamba_mf:.1f} minutes")
print(f"  Speedup: {speedup_mf:.1f}x faster")
print(f"\nAll mate-fix methods performance:")
for method, runtime in zip(sorted_methods, sorted_runtimes):
    print(f"  {method}: {runtime:.1f} minutes")