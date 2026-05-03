#!/usr/bin/env python3
"""
Convergence Analysis Plot for PN Methods with Oscillatory Corrections

This script reads the convergence data from the C++ program output and creates
a plot of log(delta_phi) vs log(eps) to verify the O(eps^(-5/2)) scaling.
"""

import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import subprocess
import sys

def run_program_and_get_data():
    """Run the C++ program and extract convergence data"""
    try:
        # Run the program and capture output
        result = subprocess.run(['./comparison_svd_TW'],
                              capture_output=True, text=True, cwd='.')

        if result.returncode != 0:
            print("Error running program:")
            print(result.stderr)
            return None

        # Extract convergence data section
        lines = result.stdout.split('\n')
        data_start = False
        data_lines = []

        for line in lines:
            if "Convergence data for log(delta_phi) vs log(eps):" in line:
                data_start = True
                continue
            elif data_start and line.strip() == "":
                continue
            elif data_start and "Example integration" in line:
                break
            elif data_start and line.strip():
                data_lines.append(line)

        return '\n'.join(data_lines[2:])  # Skip header lines

    except Exception as e:
        print(f"Error running program: {e}")
        return None

def parse_convergence_data(data_text):
    """Parse the convergence data into numpy arrays"""
    data = np.loadtxt(StringIO(data_text), usecols=(1, 2, 3, 4))
    log_eps = data[:, 0]
    log_delta_TW = data[:, 1]
    log_delta_JF = data[:, 2]
    log_delta_QLT = data[:, 3]

    return log_eps, log_delta_TW, log_delta_JF, log_delta_QLT

def create_convergence_plot(log_eps, log_delta_TW, log_delta_JF, log_delta_QLT):
    """Create the convergence plot"""
    plt.figure(figsize=(10, 8))

    # Plot the data
    plt.plot(log_eps, log_delta_TW, 'bo-', label='TW (Tucker-Will)', markersize=8, linewidth=2)
    plt.plot(log_eps, log_delta_JF, 'rs-', label='JF (Jan Fereisl)', markersize=8, linewidth=2)
    plt.plot(log_eps, log_delta_QLT, 'g^-', label='QLT (with Y-terms)', markersize=8, linewidth=2)

    # Add theoretical scaling line: slope = -5/2 = -2.5
    eps_range = np.linspace(log_eps.min(), log_eps.max(), 100)
    # Fit a line to the QLT data to show the actual scaling
    coeffs = np.polyfit(log_eps, log_delta_QLT, 1)
    slope = coeffs[0]
    intercept = coeffs[1]

    plt.plot(eps_range, slope * eps_range + intercept, 'k--',
             label='.2f', linewidth=2)

    # Add theoretical line
    theoretical_slope = -2.5  # -5/2
    # Offset to match the data roughly
    theoretical_offset = np.mean(log_delta_QLT) - theoretical_slope * np.mean(log_eps)
    plt.plot(eps_range, theoretical_slope * eps_range + theoretical_offset, 'r:',
             label='Theoretical: slope = -2.5', linewidth=2)

    # Formatting
    plt.xlabel('log(ε)', fontsize=14)
    plt.ylabel('log(δφ)', fontsize=14)
    plt.title('Convergence Analysis: Phase Error δφ vs ε\n' +
              'Evolution of p_true from 20M to 50M with oscillatory corrections\n'
              'TW & JF: Analytical methods | QLT: Numerical with Y-terms', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)

    # Add slope information
    plt.text(0.02, 0.98, f'QLT slope: {slope:.4f}\nExpected: -2.500',
             transform=plt.gca().transAxes, fontsize=11,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
    
    # Add method comparison
    tw_avg = np.mean(log_delta_TW)
    jf_avg = np.mean(log_delta_JF)
    qlt_range = log_delta_QLT.max() - log_delta_QLT.min()
    
    plt.text(0.02, 0.85, f'TW avg error: {tw_avg:.2f}\nJF avg error: {jf_avg:.2f}\nQLT range: {qlt_range:.4f}',
             transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    plt.tight_layout()
    return plt.gcf()

def save_data_to_file(log_eps, log_delta_TW, log_delta_JF, log_delta_QLT, filename='convergence_data.txt'):
    """Save the convergence data to a file for external plotting"""
    data = np.column_stack((log_eps, log_delta_TW, log_delta_JF, log_delta_QLT))
    header = "log_eps\tlog_delta_TW\tlog_delta_JF\tlog_delta_QLT"
    np.savetxt(filename, data, header=header, delimiter='\t', fmt='%.6f')
    print(f"Data saved to {filename}")

def main():
    print("Running convergence analysis...")

    # Get data from the program
    data_text = run_program_and_get_data()
    if data_text is None:
        print("Failed to get data from program")
        return

    # Parse the data
    try:
        log_eps, log_delta_TW, log_delta_JF, log_delta_QLT = parse_convergence_data(data_text)
    except Exception as e:
        print(f"Error parsing data: {e}")
        print("Raw data:")
        print(data_text)
        return

    # Create the plot
    fig = create_convergence_plot(log_eps, log_delta_TW, log_delta_JF, log_delta_QLT)

    # Save the plot
    plt.savefig('convergence_analysis.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'convergence_analysis.png'")

    # Save data to file
    save_data_to_file(log_eps, log_delta_TW, log_delta_JF, log_delta_QLT)

    # Show the plot
    plt.show()

if __name__ == "__main__":
    main()
