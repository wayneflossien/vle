#!/usr/bin/env python3
"""
Plot experimental data for Water-Glycerol binary system
"""

import numpy as np
import matplotlib.pyplot as plt

# Experimental data from the table
x1_values = [1.0000, 0.9500, 0.9000, 0.8000, 0.7000, 0.6000, 0.5000, 0.4000, 0.3000, 0.2000, 0.1000, 0.0000]
T_exp_K = [326.10, 327.10, 328.35, 331.35, 335.15, 339.90, 346.00, 353.75, 364.30, 379.65, 406.70, 497.10]
T_calc_K = [325.98, 327.07, 328.32, 331.34, 335.26, 340.35, 347.02, 355.88, 367.92, 385.19, 413.78, 496.17]

# Convert to Celsius
T_exp_C = [T - 273.15 for T in T_exp_K]
T_calc_C = [T - 273.15 for T in T_calc_K]

# Create the plot
plt.figure(figsize=(12, 8))

# Plot experimental and calculated values
plt.plot(x1_values, T_exp_C, 'ko-', linewidth=2, markersize=8, label='Experimental (Literature)')
plt.plot(x1_values, T_calc_C, 'rs-', linewidth=2, markersize=6, label='NRTL (Literature)')

plt.xlabel('Mole Fraction of Water (x₁)', fontsize=12)
plt.ylabel('Boiling Temperature (°C)', fontsize=12)
plt.title('Water-Glycerol Binary System\nBoiling Points at 13.8 kPa', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10)
plt.xlim(0, 1)

# Add text box with key information
textstr = f'Pressure: 13.8 kPa\nPure Water: {T_exp_C[0]:.1f}°C\nPure Glycerol: {T_exp_C[-1]:.1f}°C'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
         verticalalignment='top', bbox=props)

plt.tight_layout()
plt.show()

# Print the data table
print("Water-Glycerol Binary System - Experimental Data")
print("=" * 60)
print("x₁    | T_exp (°C) | T_NRTL (°C)")
print("-" * 35)

for i in range(len(x1_values)):
    print(f"{x1_values[i]:.4f} | {T_exp_C[i]:9.1f} | {T_calc_C[i]:10.1f}")

print(f"\nPure Water Boiling Point: {T_exp_C[0]:.1f}°C")
print(f"Pure Glycerol Boiling Point: {T_exp_C[-1]:.1f}°C") 