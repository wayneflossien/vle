import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from thermo import Chemical, VaporPressure, UNIFAC
from phasepy import component, mixture, vle
import warnings
warnings.filterwarnings('ignore')

class ThermoPhasePyBoilingPoint:
    """
    Calculate boiling points of binary systems using thermo and phasepy packages
    """
    
    def __init__(self, pressure=13.8):  # pressure in kPa
        self.P = pressure * 1000  # Convert to Pa
        
        # Initialize chemicals using thermo package
        self.water = Chemical('water')
        self.glycerol = Chemical('glycerol')
        
        # Create UNIFAC object for activity coefficient calculations
        self.unifac = UNIFAC.from_subgroups(['water', 'glycerol'], 
                                          [{'1': 1}, {'1': 3, '2': 3}])  # OH and CH2 groups
        
        # Initialize phasepy components using properties from thermo Chemical module
        self.water_py = component(name='water', 
                                Tc=self.water.Tc, 
                                Pc=self.water.Pc, 
                                w=self.water.omega, 
                                kappa=0.0, 
                                eos='pr')
        self.glycerol_py = component(name='glycerol', 
                                   Tc=self.glycerol.Tc, 
                                   Pc=self.glycerol.Pc, 
                                   w=self.glycerol.omega, 
                                   kappa=0.0, 
                                   eos='pr')
        
        # Create mixture for phasepy
        self.mix = mixture([self.water_py, self.glycerol_py])
        
        # Set interaction parameters (kij) - can be adjusted based on literature
        self.mix.kij = np.array([[0.0, 0.1], [0.1, 0.0]])
        
    def vapor_pressure_thermo(self, T, component_name):
        """
        Calculate vapor pressure using thermo package
        T in K, returns pressure in Pa
        """
        if component_name == 'water':
            chem = self.water
        else:
            chem = self.glycerol
            
        try:
            # Use thermo's vapor pressure calculation
            P_sat = chem.VaporPressure(T)
            return P_sat
        except:
            # Fallback to Antoine equation if thermo fails
            T_C = T - 273.15
            if component_name == 'water':
                A, B, C = 8.07131, 1730.63, 233.426
            else:  # glycerol
                A, B, C = 7.123, 1657.46, 227.02
            
            P_mmHg = 10**(A - B/(C + T_C))
            return P_mmHg * 133.322  # Convert mmHg to Pa
    
    def activity_coefficient_thermo(self, x1, T):
        """
        Calculate activity coefficients using thermo's UNIFAC
        x1: mole fraction of component 1 (water)
        """
        x2 = 1 - x1
        x = [x1, x2]
        
        try:
            # Use thermo's UNIFAC calculation
            gammas = self.unifac.gammas(x, T)
            return gammas[0], gammas[1]
        except:
            # Fallback to simplified calculation
            # For water-glycerol system, water typically shows positive deviation
            gamma1 = 1.0 + 0.5 * x2**2  # Simplified model
            gamma2 = 1.0 + 0.3 * x1**2
            return gamma1, gamma2
    
    def bubble_point_equation(self, T, x1):
        """
        Bubble point equation: P = x1*gamma1*P1_sat + x2*gamma2*P2_sat
        """
        x2 = 1 - x1
        
        # Calculate vapor pressures
        P1_sat = self.vapor_pressure_thermo(T, 'water')
        P2_sat = self.vapor_pressure_thermo(T, 'glycerol')
        
        # Calculate activity coefficients
        gamma1, gamma2 = self.activity_coefficient_thermo(x1, T)
        
        # Bubble point equation
        P_calc = x1 * gamma1 * P1_sat + x2 * gamma2 * P2_sat
        
        return P_calc - self.P
    
    def calculate_boiling_point(self, x1):
        """
        Calculate boiling point for given composition
        x1: mole fraction of water
        """
        from scipy.optimize import fsolve
        
        # Initial guess for temperature
        T_water = 273.15 + 50  # Water boiling point at 13.8 kPa ≈ 50°C
        T_glycerol = 273.15 + 200  # Glycerol boiling point at 13.8 kPa ≈ 200°C
        T_guess = x1 * T_water + (1 - x1) * T_glycerol
        
        try:
            T_bubble = fsolve(self.bubble_point_equation, T_guess, args=(x1,))[0]
            return T_bubble
        except:
            return np.nan
    
    def phasepy_vle_calculation(self, x1, T):
        """
        Use phasepy for VLE calculations
        """
        try:
            x2 = 1 - x1
            x = np.array([x1, x2])
            
            # Calculate bubble point using phasepy
            result = vle.bubble_temperature(x, self.P, self.mix, T)
            
            if result.success:
                return result.T, result.y
            else:
                return np.nan, np.nan
        except:
            return np.nan, np.nan
    
    def calculate_composition_curve(self, n_points=101):
        """
        Calculate boiling points for the entire composition range
        """
        x1_range = np.linspace(0, 1, n_points)
        boiling_points = []
        activity_coeffs = []
        vapor_compositions = []
        
        for x1 in x1_range:
            # Method 1: Using thermo package
            T_bp = self.calculate_boiling_point(x1)
            boiling_points.append(T_bp)
            
            if not np.isnan(T_bp):
                gamma1, gamma2 = self.activity_coefficient_thermo(x1, T_bp)
                activity_coeffs.append((gamma1, gamma2))
                
                # Method 2: Using phasepy for comparison
                T_phasepy, y_phasepy = self.phasepy_vle_calculation(x1, T_bp)
                if not np.isnan(T_phasepy):
                    vapor_compositions.append(y_phasepy[0])  # Water mole fraction in vapor
                else:
                    vapor_compositions.append(np.nan)
            else:
                activity_coeffs.append((np.nan, np.nan))
                vapor_compositions.append(np.nan)
        
        return x1_range, np.array(boiling_points), activity_coeffs, vapor_compositions
    
    def plot_results(self, x1_range, boiling_points, activity_coeffs=None, vapor_compositions=None):
        """
        Plot the boiling point vs composition curve and other properties
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Convert to Celsius for plotting
        T_C = boiling_points - 273.15
        
        # 1. Boiling point curve
        ax1.plot(x1_range, T_C, 'b-', linewidth=2, label='Boiling Point (Thermo)')
        ax1.set_xlabel('Mole Fraction of Water (x₁)', fontsize=12)
        ax1.set_ylabel('Boiling Point (°C)', fontsize=12)
        ax1.set_title('Water-Glycerol Binary System\nBoiling Points', fontsize=14)
        ax1.grid(True, alpha=0.3)
        
        # Add pure component boiling points
        T_water_pure = self.calculate_boiling_point(1.0) - 273.15
        T_glycerol_pure = self.calculate_boiling_point(0.0) - 273.15
        
        ax1.plot(1.0, T_water_pure, 'ro', markersize=8, label=f'Pure Water: {T_water_pure:.1f}°C')
        ax1.plot(0.0, T_glycerol_pure, 'go', markersize=8, label=f'Pure Glycerol: {T_glycerol_pure:.1f}°C')
        ax1.legend()
        
        # 2. Activity coefficients
        if activity_coeffs:
            gamma1_vals = [gamma[0] for gamma in activity_coeffs if not np.isnan(gamma[0])]
            gamma2_vals = [gamma[1] for gamma in activity_coeffs if not np.isnan(gamma[1])]
            x1_valid = [x for i, x in enumerate(x1_range) if not np.isnan(activity_coeffs[i][0])]
            
            ax2.plot(x1_valid, gamma1_vals, 'r-', linewidth=2, label='γ₁ (Water)')
            ax2.plot(x1_valid, gamma2_vals, 'g-', linewidth=2, label='γ₂ (Glycerol)')
            ax2.set_xlabel('Mole Fraction of Water (x₁)', fontsize=12)
            ax2.set_ylabel('Activity Coefficient (γ)', fontsize=12)
            ax2.set_title('Activity Coefficients (UNIFAC)', fontsize=14)
            ax2.grid(True, alpha=0.3)
            ax2.legend()
        
        # 3. Vapor composition (y-x diagram)
        if vapor_compositions:
            y_valid = [y for y in vapor_compositions if not np.isnan(y)]
            x1_valid = [x for i, x in enumerate(x1_range) if not np.isnan(vapor_compositions[i])]
            
            if len(y_valid) > 0:
                ax3.plot(x1_valid, y_valid, 'm-', linewidth=2, label='Vapor Composition')
                ax3.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='y = x')
                ax3.set_xlabel('Liquid Mole Fraction of Water (x₁)', fontsize=12)
                ax3.set_ylabel('Vapor Mole Fraction of Water (y₁)', fontsize=12)
                ax3.set_title('Vapor-Liquid Equilibrium (y-x)', fontsize=14)
                ax3.grid(True, alpha=0.3)
                ax3.legend()
        
        # 4. Vapor pressure comparison
        T_range = np.linspace(273.15 + 20, 273.15 + 150, 50)
        P_water = [self.vapor_pressure_thermo(T, 'water')/1000 for T in T_range]
        P_glycerol = [self.vapor_pressure_thermo(T, 'glycerol')/1000 for T in T_range]
        
        ax4.plot(T_range - 273.15, P_water, 'b-', linewidth=2, label='Water')
        ax4.plot(T_range - 273.15, P_glycerol, 'g-', linewidth=2, label='Glycerol')
        ax4.axhline(y=self.P/1000, color='r', linestyle='--', alpha=0.7, label=f'P = {self.P/1000:.1f} kPa')
        ax4.set_xlabel('Temperature (°C)', fontsize=12)
        ax4.set_ylabel('Vapor Pressure (kPa)', fontsize=12)
        ax4.set_title('Vapor Pressure vs Temperature', fontsize=14)
        ax4.set_yscale('log')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        
        plt.tight_layout()
        plt.show()
    
    def save_results(self, x1_range, boiling_points, activity_coeffs=None, 
                    vapor_compositions=None, filename='thermo_phasepy_results.csv'):
        """
        Save results to CSV file
        """
        data = {
            'Mole_Fraction_Water': x1_range,
            'Mole_Fraction_Glycerol': 1 - x1_range,
            'Boiling_Point_K': boiling_points,
            'Boiling_Point_C': boiling_points - 273.15
        }
        
        if activity_coeffs:
            gamma1_vals = [gamma[0] if not np.isnan(gamma[0]) else np.nan for gamma in activity_coeffs]
            gamma2_vals = [gamma[1] if not np.isnan(gamma[1]) else np.nan for gamma in activity_coeffs]
            data['Activity_Coefficient_Water'] = gamma1_vals
            data['Activity_Coefficient_Glycerol'] = gamma2_vals
        
        if vapor_compositions:
            y_vals = [y if not np.isnan(y) else np.nan for y in vapor_compositions]
            data['Vapor_Mole_Fraction_Water'] = y_vals
        
        results_df = pd.DataFrame(data)
        results_df.to_csv(filename, index=False)
        print(f"Results saved to {filename}")
        return results_df
    
    def compare_methods(self, x1=0.5):
        """
        Compare different calculation methods
        """
        print(f"Comparison of Methods for x₁ = {x1:.2f}")
        print("=" * 50)
        
        # Method 1: Thermo package
        T_thermo = self.calculate_boiling_point(x1)
        gamma1_thermo, gamma2_thermo = self.activity_coefficient_thermo(x1, T_thermo)
        
        print(f"Thermo Package:")
        print(f"  Boiling Point: {T_thermo - 273.15:.1f}°C")
        print(f"  Activity Coefficients: γ₁ = {gamma1_thermo:.3f}, γ₂ = {gamma2_thermo:.3f}")
        
        # Method 2: Phasepy
        T_phasepy, y_phasepy = self.phasepy_vle_calculation(x1, T_thermo)
        
        if not np.isnan(T_phasepy):
            print(f"\nPhasepy Package:")
            print(f"  Bubble Temperature: {T_phasepy - 273.15:.1f}°C")
            print(f"  Vapor Composition: y₁ = {y_phasepy[0]:.3f}, y₂ = {y_phasepy[1]:.3f}")
        else:
            print(f"\nPhasepy Package: Calculation failed")
        
        # Method 3: Pure component properties
        print(f"\nPure Component Properties (from Thermo):")
        print(f"  Water Critical Temperature: {self.water.Tc - 273.15:.1f}°C")
        print(f"  Water Critical Pressure: {self.water.Pc/1e6:.1f} MPa")
        print(f"  Water Acentric Factor: {self.water.omega:.3f}")
        print(f"  Glycerol Critical Temperature: {self.glycerol.Tc - 273.15:.1f}°C")
        print(f"  Glycerol Critical Pressure: {self.glycerol.Pc/1e6:.1f} MPa")
        print(f"  Glycerol Acentric Factor: {self.glycerol.omega:.3f}")
        
        # Display phasepy component properties
        print(f"\nPhasepy Component Properties (from Thermo):")
        print(f"  Water (Phasepy): Tc={self.water_py.Tc:.1f}K, Pc={self.water_py.Pc/1e6:.1f}MPa, ω={self.water_py.w:.3f}")
        print(f"  Glycerol (Phasepy): Tc={self.glycerol_py.Tc:.1f}K, Pc={self.glycerol_py.Pc/1e6:.1f}MPa, ω={self.glycerol_py.w:.3f}")
    
    def display_component_properties(self):
        """
        Display comprehensive component properties from thermo Chemical module
        """
        print("Component Properties from Thermo Chemical Module")
        print("=" * 60)
        
        # Water properties
        print(f"\nWater Properties:")
        print(f"  Critical Temperature: {self.water.Tc:.2f} K ({self.water.Tc - 273.15:.1f}°C)")
        print(f"  Critical Pressure: {self.water.Pc/1e6:.2f} MPa ({self.water.Pc/1e5:.1f} bar)")
        print(f"  Critical Volume: {self.water.Vc*1000:.2f} L/mol")
        print(f"  Acentric Factor: {self.water.omega:.3f}")
        print(f"  Molecular Weight: {self.water.MW:.2f} g/mol")
        
        # Glycerol properties
        print(f"\nGlycerol Properties:")
        print(f"  Critical Temperature: {self.glycerol.Tc:.2f} K ({self.glycerol.Tc - 273.15:.1f}°C)")
        print(f"  Critical Pressure: {self.glycerol.Pc/1e6:.2f} MPa ({self.glycerol.Pc/1e5:.1f} bar)")
        print(f"  Critical Volume: {self.glycerol.Vc*1000:.2f} L/mol")
        print(f"  Acentric Factor: {self.glycerol.omega:.3f}")
        print(f"  Molecular Weight: {self.glycerol.MW:.2f} g/mol")
        
        # Phasepy component properties
        print(f"\nPhasepy Component Properties:")
        print(f"  Water: Tc={self.water_py.Tc:.1f}K, Pc={self.water_py.Pc/1e6:.1f}MPa, ω={self.water_py.w:.3f}")
        print(f"  Glycerol: Tc={self.glycerol_py.Tc:.1f}K, Pc={self.glycerol_py.Pc/1e6:.1f}MPa, ω={self.glycerol_py.w:.3f}")
        
        # Verify that phasepy properties match thermo properties
        print(f"\nProperty Verification:")
        print(f"  Water Tc match: {abs(self.water.Tc - self.water_py.Tc) < 0.1}")
        print(f"  Water Pc match: {abs(self.water.Pc - self.water_py.Pc) < 1000}")
        print(f"  Water ω match: {abs(self.water.omega - self.water_py.w) < 0.001}")
        print(f"  Glycerol Tc match: {abs(self.glycerol.Tc - self.glycerol_py.Tc) < 0.1}")
        print(f"  Glycerol Pc match: {abs(self.glycerol.Pc - self.glycerol_py.Pc) < 1000}")
        print(f"  Glycerol ω match: {abs(self.glycerol.omega - self.glycerol_py.w) < 0.001}")

def main():
    """
    Main function to run the thermo/phasepy boiling point calculations
    """
    print("Water-Glycerol Binary System Boiling Point Calculator")
    print("Using Thermo and Phasepy Packages at 13.8 kPa")
    print("=" * 60)
    
    # Initialize calculator
    calculator = ThermoPhasePyBoilingPoint(pressure=13.8)
    
    # Display component properties from thermo Chemical module
    calculator.display_component_properties()
    
    # Compare methods for a specific composition
    calculator.compare_methods(x1=0.5)
    
    # Calculate boiling points for different compositions
    print("\nCalculating boiling points and activity coefficients...")
    x1_range, boiling_points, activity_coeffs, vapor_compositions = calculator.calculate_composition_curve()
    
    # Display key results
    print("\nKey Results:")
    print(f"Pure Water Boiling Point: {boiling_points[-1] - 273.15:.1f}°C")
    print(f"Pure Glycerol Boiling Point: {boiling_points[0] - 273.15:.1f}°C")
    
    # Find minimum boiling point
    valid_indices = ~np.isnan(boiling_points)
    if np.any(valid_indices):
        min_idx = np.argmin(boiling_points[valid_indices])
        actual_min_idx = np.where(valid_indices)[0][min_idx]
        print(f"Minimum Boiling Point: {boiling_points[actual_min_idx] - 273.15:.1f}°C at x₁ = {x1_range[actual_min_idx]:.3f}")
    
    # Plot results
    calculator.plot_results(x1_range, boiling_points, activity_coeffs, vapor_compositions)
    
    # Save results
    results_df = calculator.save_results(x1_range, boiling_points, activity_coeffs, vapor_compositions)
    
    # Display sample results
    print("\nSample Results:")
    print(results_df.iloc[::20])  # Show every 20th result
    
    return results_df

if __name__ == "__main__":
    results = main() 