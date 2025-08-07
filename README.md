# Water-Glycerol Binary System Boiling Point Calculator

This project provides a comprehensive tool for calculating boiling points of water-glycerol binary mixtures using the UNIFAC (UNIQUAC Functional-group Activity Coefficients) method at different compositions and pressures.

## Features

- **UNIFAC Method Implementation**: Uses the UNIFAC group contribution method for accurate activity coefficient calculations
- **Antoine Equation**: Implements vapor pressure calculations using Antoine equation parameters
- **Composition Range**: Calculates boiling points from 0 to 1 mole fraction of water
- **Pressure Flexibility**: Configurable pressure (default: 13.8 kPa)
- **Visualization**: Generates plots of boiling points vs composition and activity coefficients
- **Data Export**: Saves results to CSV files for further analysis
- **Azeotrope Detection**: Identifies minimum boiling azeotropes if present
- **Experimental Comparison**: Compares calculated values with experimental data from literature

## Files Description

1. **`thermo_phasepy_boiling_point.py`**: Advanced implementation using thermo and phasepy packages for accurate thermodynamic calculations
2. **`plot_experimental_data.py`**: Simple script to plot experimental data from literature
3. **`requirements.txt`**: Python dependencies
4. **`README.md`**: This documentation file

## Installation

1. Clone or download the project files
2. Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

#### Advanced Thermo/Phasepy Implementation (Recommended)
Run the advanced calculation script:

```bash
python thermo_phasepy_boiling_point.py
```

This will:
- Calculate boiling points for the entire composition range (0-1 mole fraction water)
- Generate comprehensive plots including boiling points, activity coefficients, vapor compositions, and vapor pressures
- Save results to `thermo_phasepy_results.csv`
- Display key results including pure component boiling points
- Compare different calculation methods (thermo vs phasepy)

#### Experimental Data Plot
Run the experimental data plot script:

```bash
python plot_experimental_data.py
```

This will:
- Plot experimental boiling points from literature
- Show NRTL calculated values
- Display the exact data points from the experimental table
- Print a detailed comparison table

### Custom Usage

#### Advanced Thermo/Phasepy Implementation
```python
from thermo_phasepy_boiling_point import ThermoPhasePyBoilingPoint

# Initialize calculator at 13.8 kPa
calculator = ThermoPhasePyBoilingPoint(pressure=13.8)

# Calculate boiling point for 50% water composition
x_water = 0.5
T_bp = calculator.calculate_boiling_point(x_water)
T_bp_C = T_bp - 273.15
print(f"Boiling point at {x_water:.1%} water: {T_bp_C:.1f}°C")

# Calculate activity coefficients using thermo package
gamma1, gamma2 = calculator.activity_coefficient_thermo(x_water, T_bp)
print(f"Activity coefficients: γ₁ = {gamma1:.3f}, γ₂ = {gamma2:.3f}")

# Use phasepy for VLE calculations
T_phasepy, y_phasepy = calculator.phasepy_vle_calculation(x_water, T_bp)
print(f"Vapor composition: y₁ = {y_phasepy[0]:.3f}, y₂ = {y_phasepy[1]:.3f}")

# Compare different methods
calculator.compare_methods(x1=x_water)
```

## Method Details

### UNIFAC Method

The UNIFAC method calculates activity coefficients using group contributions:

1. **Combinatorial Part**: Accounts for molecular size and shape differences
2. **Residual Part**: Accounts for group-group interactions

### Group Parameters

- **Water**: 1 OH group
- **Glycerol**: 3 OH groups, 3 CH2 groups

### Antoine Equation

Vapor pressure calculation:
```
log₁₀(P) = A - B/(C + T)
```
Where:
- P = vapor pressure (mmHg)
- T = temperature (°C)
- A, B, C = component-specific parameters

### Bubble Point Equation

```
P = x₁·γ₁·P₁^sat + x₂·γ₂·P₂^sat
```
Where:
- P = total pressure
- xᵢ = mole fraction of component i
- γᵢ = activity coefficient of component i
- Pᵢ^sat = saturation vapor pressure of component i

## Output Files

### CSV Results

The calculator generates CSV files with columns:
- `Mole_Fraction_Water`: Mole fraction of water (0-1)
- `Mole_Fraction_Glycerol`: Mole fraction of glycerol (0-1)
- `Boiling_Point_K`: Boiling point in Kelvin
- `Boiling_Point_C`: Boiling point in Celsius
- `Activity_Coefficient_Water`: Activity coefficient of water
- `Activity_Coefficient_Glycerol`: Activity coefficient of glycerol
- `Vapor_Mole_Fraction_Water`: Vapor composition of water (if available)

### Plots

1. **Boiling Point vs Composition**: Shows how boiling point varies with water concentration
2. **Activity Coefficients**: Shows activity coefficients vs composition
3. **Vapor-Liquid Equilibrium**: Shows vapor composition vs liquid composition
4. **Vapor Pressure**: Shows vapor pressure vs temperature for pure components

## Key Results

At 13.8 kPa pressure:
- **Pure Water Boiling Point**: ~53°C (experimental)
- **Pure Glycerol Boiling Point**: ~224°C (experimental)
- **Minimum Boiling Point**: Varies with composition (azeotrope analysis)

## Experimental Data Comparison

The project includes experimental data from literature for validation:
- **Data Points**: 12 compositions from 0 to 1 mole fraction water
- **Temperature Range**: 53°C to 224°C
- **NRTL Model Accuracy**: Average deviation < 0.5%
- **Maximum Deviation**: 1.74% at low water concentrations

## Technical Notes

### Assumptions

1. **Ideal Gas Phase**: Vapor phase is assumed to be ideal
2. **Liquid Phase Non-ideality**: Accounted for using UNIFAC activity coefficients
3. **Temperature Range**: Antoine equation parameters are valid for specific temperature ranges
4. **Group Contributions**: Simplified group interaction parameters are used

### Limitations

1. **Temperature Range**: Antoine equation may not be accurate outside its parameter range
2. **Group Parameters**: Simplified UNIFAC parameters may not capture all molecular interactions
3. **Pressure Effects**: High pressure effects on activity coefficients are not included

### Accuracy

The method provides reasonable estimates for:
- Boiling point trends vs composition
- Activity coefficient behavior
- Azeotrope identification

For high-precision applications, experimental validation is recommended.

## Dependencies

### Advanced Implementation (Thermo/Phasepy)
- **thermo**: Comprehensive thermodynamic property calculations
- **phasepy**: Phase equilibrium calculations using cubic equations of state
- **numpy**: Numerical computations
- **matplotlib**: Plotting and visualization
- **scipy**: Optimization (fsolve for bubble point calculation)
- **pandas**: Data handling and CSV export

## References

1. Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975). Group-contribution estimation of activity coefficients in nonideal liquid mixtures. AIChE Journal, 21(6), 1086-1099.
2. Gmehling, J., Li, J., & Schiller, M. (1993). A modified UNIFAC model. 2. Present parameter matrix and results for different thermodynamic properties. Industrial & Engineering Chemistry Research, 32(1), 178-193.
3. Antoine, C. (1888). Tensions des vapeurs; nouvelle relation entre les tensions et les températures. Comptes Rendus des Séances de l'Académie des Sciences, 107, 778-780.
4. Bell, I. H., Wronski, J., Quoilin, S., & Lemort, V. (2014). Pure and pseudo-pure fluid thermophysical property evaluation and the open-source thermophysical property library CoolProp. Industrial & Engineering Chemistry Research, 53(6), 2498-2508.
5. Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001). The properties of gases and liquids. McGraw-Hill Education.

## License

This code is provided for educational and research purposes. Please cite appropriate references when using this work in publications. 