#!/usr/bin/env python3
"""
Plot showing percentage improvement of JF over TW as a function of epsilon
For specific epsilon values: 0.1, 0.01, 0.001, 0.0001 to 0.00000001
"""

import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np

def parse_coefficient_section(output):
    """Parse the PN coefficient comparison section to extract TW vs JF relDiff for each order"""
    tw_values = []
    jf_values = []
    orders = []
    
    lines = output.split('\n')
    
    # Find the start of the dp/dtheta section
    start_idx = None
    for i, line in enumerate(lines):
        if ' dp/dtheta:' in line:
            start_idx = i
            break
    
    if start_idx is None:
        return None
    
    # Extract the data lines (they contain PN orders)
    for i in range(start_idx + 1, len(lines)):
        line = lines[i]
        
        # Stop at the de/dtheta line
        if 'de/dtheta:' in line:
            break
        
        # Look for lines with PN orders
        if '2.5PN' in line or '3.5PN' in line or '4PN' in line or '4.5PN' in line:
            # Split and remove empty strings
            parts = [p for p in line.split() if p]
            
            # The structure is: PNorder, value, ±, uncertainty, TW_val, JF_val, relDiff_TW, relDiff_JF, Winner, *marker
            # So: parts[0]=order, parts[6]=relDiff(QLT,TW), parts[7]=relDiff(QLT,JF)
            if len(parts) >= 8:
                try:
                    order = parts[0].strip()
                    rd_tw = float(parts[6])
                    rd_jf = float(parts[7])
                    orders.append(order)
                    tw_values.append(rd_tw)
                    jf_values.append(rd_jf)
                except (ValueError, IndexError):
                    pass
    
    if not tw_values or not jf_values:
        return None
    
    return {
        'orders': orders,
        'tw_values': tw_values,
        'jf_values': jf_values
    }

def parse_total_reldiff(output):
    """Parse total relative differences from Section 2"""
    # Look for pattern like "Total relDiff\s+TW=... JF=..."
    patterns = [
        r'Total relDiff\s+TW=([\d.eE+-]+)\s+JF=([\d.eE+-]+)',
        r'TOTAL.*?([\d.eE+-]+).*?([\d.eE+-]+).*?([\d.eE+-]+)',
    ]
    
    for pattern in patterns:
        match = re.search(pattern, output)
        if match:
            try:
                tw_val = float(match.group(1))
                jf_val = float(match.group(2))
                return tw_val, jf_val
            except (ValueError, IndexError):
                pass
    
    return None, None

# Run both programs
print("Running comparison_svd_TW...")
result_svd = subprocess.run(['./comparison_svd_TW'], capture_output=True, text=True, cwd='.')
output_svd = result_svd.stdout

print("\nParsing outputs...")

# Parse coefficient data (shows TW vs JF at each PN order)
svd_coeff = parse_coefficient_section(output_svd)

# Parse total relative differences
svd_tw_total, svd_jf_total = parse_total_reldiff(output_svd)

print(f"SVD - TW: {svd_tw_total}, JF: {svd_jf_total}")

# Create visualization
fig, ax = plt.subplots(figsize=(12, 6))
fig.suptitle('SVD Method: TW vs JF Performance by PN Order', fontsize=14, fontweight='bold')

# Plot: Relative differences at each PN order
if svd_coeff:
    x_pos = np.arange(len(svd_coeff['orders']))
    width = 0.35
    bars1 = ax.bar(x_pos - width/2, svd_coeff['tw_values'], width, label='TW RelDiff', alpha=0.8, color='steelblue', edgecolor='black', linewidth=1)
    bars2 = ax.bar(x_pos + width/2, svd_coeff['jf_values'], width, label='JF RelDiff', alpha=0.8, color='coral', edgecolor='black', linewidth=1)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0.001:  # Only label significant values
                ax.text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.4f}',
                        ha='center', va='bottom', fontsize=9)
    
    ax.set_xlabel('PN Order', fontsize=12, fontweight='bold')
    ax.set_ylabel('Relative Difference', fontsize=12, fontweight='bold')
    ax.set_title('SVD Method: dp/dtheta Coefficient Comparison', fontsize=12, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(svd_coeff['orders'], fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_yscale('log')
else:
    ax.text(0.5, 0.5, 'No data parsed for SVD', ha='center', va='center', transform=ax.transAxes)

plt.tight_layout()
plt.savefig('jf_vs_tw_comparison.png', dpi=150, bbox_inches='tight')
print("✓ Saved: jf_vs_tw_comparison.png")

# Calculate improvement percentages
print("\n" + "="*80)
print("JF IMPROVEMENT OVER TW")
print("="*80)

if svd_tw_total and svd_jf_total:
    svd_improvement = ((svd_tw_total - svd_jf_total) / abs(svd_tw_total)) * 100
    svd_better = "Yes" if svd_jf_total < svd_tw_total else "No"
    print(f"\nSVD METHOD:")
    print(f"  TW total relative difference: {svd_tw_total:.16f}")
    print(f"  JF total relative difference: {svd_jf_total:.16f}")
    print(f"  JF is better: {svd_better}")
    print(f"  Improvement: {svd_improvement:.16f}%")

print("="*80)

#plt.show()
