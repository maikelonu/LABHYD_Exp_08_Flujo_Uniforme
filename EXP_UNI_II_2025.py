import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from statsmodels.formula.api import ols

# Introduction
print("\n//////////////////////////////////////////////////////////////////////////////////")
print("INSTITUTO TECNOLOGICO DE COSTA RICA")
print("Escuela de Ingeniería en Construcción")
print("https://www.tec.ac.cr")
print("Session: FLUJO UNIFORME")

print("\nM.Sc. Eng. Maikel Méndez M")
print("Water Resources + GIS + DataScience")
print("Instituto Tecnológico de Costa Rica")
print("https://www.tec.ac.cr")
print("https://orcid.org/0000-0003-1919-141X")
print("https://www.scopus.com/authid/detail.uri?authorId=51665581300")
print("https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en")
print("https://www.youtube.com/c/maikelmendez")
print("https://github.com/maikelonu")
print("//////////////////////////////////////////////////////////////////////////////////")

print("\n# INFO:")
print("Análisis gráfico avanzado")
print("Normalización y homogenización de variables")
print("Exportación ASCII")
print("//////////////////////////////////////////////////////////////////////////////////")

# Scientific notation is suppress
print("\n# Suppress scientific notation for clarity")
pd.set_option('display.float_format', lambda x: '%.3f' % x)

# Set the working directory
os.chdir("/home/shoe/Dropbox/Academics/LAB_Esencial/PYTHON_LABHYD_Exp/PYTHON_LABHYD_Exp_08")
print("\n# Working directory is selected:")
print(os.getcwd())

print("\n/////////////////////////////////////////////////////////////")
print("BLOCK: Custom function, round_df to specific digits")
print("/////////////////////////////////////////////////////////////")
def round_df(df, digits):
    return df.round(digits)

print("\n////////////////////////////////////////////////////////")
print("BLOCK: Declarations")
print("////////////////////////////////////////////////////////")
base_m = 0.086
print("base_m <- 0.086")
visco = 1e-06
print("visco <- 1e-06")

print("\n////////////////////////////////////////////////////////")
print("BLOCK: Data input")
print("////////////////////////////////////////////////////////")

print("\n# Input data is loaded and a data.frame is created")
df_base = pd.read_csv("uniforme.txt", sep="\t")

# Function to calculate detailed descriptive statistics for each column in the DataFrame
def detailed_desc(df):
    desc = {}  # Create an empty dictionary to store statistics for each column
    for col in df.columns:  # Loop through each column in the DataFrame
        if pd.api.types.is_numeric_dtype(df[col]):  # Check if the column is numeric
            n = len(df[col])  # Get the number of observations in the column
            mean = df[col].mean()  # Calculate the mean of the column
            std_dev = df[col].std()  # Calculate the standard deviation of the column
            vcoef = std_dev / mean if mean != 0 else np.nan  # Calculate the coefficient of variation (standard deviation divided by mean)
            mad = np.mean(np.abs(df[col] - mean))  # Calculate the mean absolute deviation (average distance from the mean)
            q1 = df[col].quantile(0.25)  # Calculate the 1st quartile (25th percentile)
            q3 = df[col].quantile(0.75)  # Calculate the 3rd quartile (75th percentile)
            iqr = q3 - q1  # Calculate the interquartile range (difference between 3rd and 1st quartiles)
            median = df[col].median()  # Calculate the median (middle value) of the column
            range_val = df[col].max() - df[col].min()  # Calculate the range (difference between max and min values)
            skew = df[col].skew()  # Calculate the skewness (measure of asymmetry) of the column
            kurt = df[col].kurt()  # Calculate the kurtosis (measure of "tailedness") of the column
            # Calculate the lower and upper bounds of the 95% confidence interval for the mean
            mean_ci_lower = mean - 1.96 * std_dev / np.sqrt(n)
            mean_ci_upper = mean + 1.96 * std_dev / np.sqrt(n)
            # Store all the calculated statistics in the dictionary under the column name
            desc[col] = {
                'length': n,
                'mean': mean,
                'std_dev': std_dev,
                'vcoef': vcoef,
                'mad': mad,
                'q1': q1,
                'median': median,
                'q3': q3,
                'range': range_val,
                'skew': skew,
                'kurt': kurt,
                'mean_ci_lower': mean_ci_lower,
                'mean_ci_upper': mean_ci_upper,
                'value_freq': df[col].value_counts().sort_index().to_dict()  # Frequency of each value in the column
            }
        else:  # If the column is not numeric (e.g., categorical)
            # Store summary statistics for non-numeric columns
            desc[col] = {
                'length': len(df[col]),  # Number of observations
                'unique': df[col].nunique(),  # Number of unique values
                'levels': df[col].unique(),  # List of unique values
                'dupes': df[col].duplicated().any()  # Check if there are any duplicates
            }
    return desc  # Return the dictionary containing the descriptive statistics

# Generate the detailed descriptive statistics
desc_stats = detailed_desc(df_base)

# Function to print the detailed descriptive statistics
def print_detailed_desc(desc_stats):
    print("\n--------------------------------------------------------------------------------\n")
    print(f"Describe df.base (data.frame):\n")
    print(f"data frame:\t{len(df_base)} obs. of {df_base.shape[1]} variables")
    print(f"\t\t{len(df_base.dropna())} complete cases ({(len(df_base.dropna()) / len(df_base)) * 100:.1f}%)")
    print("\n  Nr  ColName     Class        NAs  Levels")
    
    for i, (col, stats) in enumerate(desc_stats.items(), 1):
        col_class = str(df_base[col].dtype)  # Convert the dtype to a string
        print(f"  {i}   {col:<10}  {col_class:<10}    .")
    
    for i, (col, stats) in enumerate(desc_stats.items(), 1):
        print("\n--------------------------------------------------------------------------------\n")
        print(f"{i} - {col} ({str(df_base[col].dtype)})\n")  # Convert the dtype to a string here as well
        if pd.api.types.is_numeric_dtype(df_base[col]):
            print(f"  length       n    NAs  unique    0s  mean  meanCI'")
            print(f"      {stats['length']}      {stats['length']}      0       {df_base[col].nunique()}     0  {stats['mean']:.2f}    {stats['mean_ci_lower']:.2f}")
            print(f"          100.0%   0.0%          0.0%          {stats['mean_ci_upper']:.2f}")
            print("                                                    ")
            print(f"     .05     .10    .25  median   .75   .90     .95")
            print(f"    {df_base[col].quantile(0.05):.2f}    {df_base[col].quantile(0.10):.2f}   {stats['q1']:.2f}    {stats['median']:.2f}  {stats['q3']:.2f}  {df_base[col].quantile(0.90):.2f}    {df_base[col].quantile(0.95):.2f}")
            print("                                                    ")
            print(f"   range      sd  vcoef     mad   IQR  skew    kurt")
            print(f"   {stats['range']:.2f}    {stats['std_dev']:.2f}   {stats['vcoef']:.2f}    {stats['mad']:.2f}  {stats['q3'] - stats['q1']:.2f}  {stats['skew']:.2f}    {stats['kurt']:.2f}")
            print("                                                    ")
            print("   value  freq   perc  cumfreq  cumperc")
            cumfreq = 0
            cumperc = 0
            for j, (val, freq) in enumerate(stats['value_freq'].items(), 1):
                perc = freq / stats['length'] * 100
                cumfreq += freq
                cumperc += perc
                print(f"{j:2d}    {val:<8}  {freq:5}  {perc:6.1f}%  {cumfreq:8}  {cumperc:8.1f}%")
        else:
            print(f"  length      n    NAs unique levels  dupes")
            print(f"      {stats['length']}      {stats['length']}      0      {stats['unique']}      {len(stats['levels'])}      {'y' if stats['dupes'] else 'n'}")
            print("                                                    ")
            print("   level  freq   perc  cumfreq  cumperc")
            freq_table = df_base[col].value_counts().to_dict()
            cumfreq = 0
            cumperc = 0
            for j, (level, freq) in enumerate(freq_table.items(), 1):
                perc = freq / stats['length'] * 100
                cumfreq += freq
                cumperc += perc
                print(f"{j:2d}    {level:<6}  {freq:5}  {perc:6.1f}%  {cumfreq:8}  {cumperc:8.1f}%")

# Print the detailed descriptive statistics
print_detailed_desc(desc_stats)

# Function to retrieve and print column names in the desired format
def print_column_names(df):
    print("\n# names {base} function is requested")
    print("names(df.base)")
    column_names = df.columns.tolist()
    formatted_column_names = " ".join(f'"{name}"' for name in column_names)
    print(f"[1] {formatted_column_names}\n")

# Call the function to print the column names
print_column_names(df_base)

# Variables
# A new slope factor variable (as %) is created
df_base['slope_perc'] = (df_base['slope_m_m'] * 100).astype(str)
# Flow units are converted from m3/h to m3/s
df_base['q_m3_s'] = df_base['q_m3_h'] / 3600
# Water depth units are converted from cm to m
df_base['y_m'] = df_base['y_cm'] / 100
# Hydraulic area is calculated
df_base['area'] = df_base['y_m'] * base_m
# Hydraulic perimeter is calculated
df_base['perimeter'] = (df_base['y_m'] * 2) + base_m
# Hydraulic radius is calculated
df_base['radius'] = df_base['area'] / df_base['perimeter']
# Square root of hydraulic radius is calculated
df_base['radius.root'] = np.sqrt(df_base['radius'])
# Water velocity is calculated
df_base['vel'] = df_base['q_m3_s'] / df_base['area']
# Chezy experimental coefficient is calculated
df_base['c.exp'] = df_base['vel'] / np.sqrt(df_base['radius'] * df_base['slope_m_m'])
# n of Manning is calculated
df_base['n.Manning'] = (df_base['radius'] ** (1/6)) / df_base['c.exp']
# m of Kutter is calculated
df_base['m.Kutter'] = (df_base['radius.root'] * 100 / df_base['c.exp']) - df_base['radius.root']
# s of Bazin is calculated
df_base['s.Bazin'] = ((87 / df_base['c.exp']) - 1) * df_base['radius.root']
# Froude number is calculated
df_base['Froude'] = df_base['vel'] / np.sqrt(df_base['area'] * 9.81 / base_m)
# Coefficients are normalized based on their mean
df_base['n.Manning.norm'] = df_base['n.Manning'] / df_base['n.Manning'].median()
df_base['m.Kutter.norm'] = df_base['m.Kutter'] / df_base['m.Kutter'].median()
df_base['s.Bazin.norm'] = df_base['s.Bazin'] / df_base['s.Bazin'].median()

# Plots
# Function to create combined plots for a given variable in the DataFrame
def create_combined_plot(df, variable, color, plot_number, dtype):
    # Set up a grid of subplots with different height ratios
    fig, axs = plt.subplots(3, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1, 1]})

    # Plot the histogram with KDE (Kernel Density Estimation)
    sns.histplot(df[variable], kde=True, color=color, ax=axs[0])
    axs[0].set_title(f"{plot_number}. {variable} ({dtype})", fontsize=18)
    axs[0].set_xlabel("")
    axs[0].set_ylabel("")

    # Plot the horizontal boxplot to show spread and outliers
    sns.boxplot(x=df[variable], ax=axs[1], color='lightgray')
    axs[1].set_xlabel("")
    axs[1].set_ylabel("")

    # Plot the empirical CDF (Cumulative Distribution Function) as a step plot
    sorted_values = np.sort(df[variable])
    step_y = np.arange(1, len(sorted_values) + 1) / len(sorted_values)
    axs[2].step(sorted_values, step_y, where='mid', color=color)
    axs[2].set_ylim(0, 1)
    axs[2].set_yticks(np.arange(0, 1.25, 0.25))
    axs[2].set_xlabel("")
    axs[2].set_ylabel("")

    # Adjust layout to prevent overlapping labels
    plt.tight_layout()
    plt.show()
    print(f"Plot {plot_number} for {variable} ({dtype}) has been generated.")

# Generate combined plots for the four columns with respective colors and data types
create_combined_plot(df_base, 'q_m3_h', 'blue', 1, 'int64')
create_combined_plot(df_base, 'y_cm', 'green', 2, 'float64')
create_combined_plot(df_base, 'slope_m_m', 'red', 3, 'float64')
create_combined_plot(df_base, 'group', 'purple', 4, 'category')

# Set up the plotting style
sns.set_theme(style="whitegrid", context="talk")

# Round the velocity values to 1 decimal place for display
df_base['tempV'] = round(df_base['vel'], 2)

# Round the relevant DataFrames to 3 decimal places
df_output = df_base.round(3)

# Save the output DataFrames to CSV files
df_output.to_csv("df.output.csv", index=False)

# Plot 1: Boxplot of water depth (y_cm) by slope percentage (slope_perc)
# This plot includes labels for flow rate (q_m3_h) and velocity (tempV)
plt.figure(figsize=(10, 8))
sns.boxplot(x='slope_perc', y='y_cm', data=df_base, width=0.75, palette='Set3')
sns.scatterplot(x='slope_perc', y='y_cm', data=df_base, s=100, legend=False)

# Add text labels for flow rate and velocity
for i in range(df_base.shape[0]):
    plt.text(x=i % len(df_base['slope_perc'].unique()), y=df_base['y_cm'][i], s=df_base['q_m3_h'][i],
             color='black', ha='right', va='bottom', fontsize=12, rotation=45)
    plt.text(x=i % len(df_base['slope_perc'].unique()), y=df_base['y_cm'][i], s=df_base['tempV'][i],
             color='blue', ha='right', va='top', fontsize=10)

# Set the title and labels for the axes
plt.title("Boxplot - Distribution of Water Depth by Slope. Labeled by Flow Rate and Velocity", fontsize=16)
plt.xlabel("Slope (%)", fontsize=14)
plt.ylabel("Water Depth (cm)", fontsize=14)
plt.xticks(rotation=0)  # Ensure x-axis labels are horizontal
plt.grid(True)  # Add grid lines
plt.show()  # Display the plot
print("\nPlot 1: Boxplot of Water Depth by Slope with Flow Rate and Velocity labels has been generated.")

# Remove the temporary velocity column
df_base.drop(columns=['tempV'], inplace=True)

# Plot 2: Boxplot of water depth (y_cm) by slope percentage (slope_perc)
# Labeled by Froude number (tempF) and velocity (tempV)
df_base['tempF'] = round(df_base['Froude'], 2)  # Round Froude number to 1 decimal place
df_base['tempV'] = round(df_base['vel'], 2)  # Round velocity to 1 decimal place

plt.figure(figsize=(10, 8))
sns.boxplot(x='slope_perc', y='y_cm', data=df_base, width=0.75, palette='Set3')
sns.scatterplot(x='slope_perc', y='y_cm', data=df_base, s=100, legend=False)

# Add text labels for Froude number and velocity
for i in range(df_base.shape[0]):
    plt.text(x=i % len(df_base['slope_perc'].unique()), y=df_base['y_cm'][i], s=df_base['tempF'][i],
             color='black', ha='right', va='bottom', fontsize=12, rotation=45)
    plt.text(x=i % len(df_base['slope_perc'].unique()), y=df_base['y_cm'][i], s=df_base['tempV'][i],
             color='blue', ha='right', va='top', fontsize=10)

# Set the title and labels for the axes
plt.title("Boxplot - Distribution of Water Depth by Slope. Labeled by Froude Number and Velocity", fontsize=16)
plt.xlabel("Slope (%)", fontsize=14)
plt.ylabel("Water Depth (cm)", fontsize=14)
plt.xticks(rotation=0)
plt.grid(True)
plt.show()
print("Plot 2: Boxplot of Water Depth by Slope with Froude Number and Velocity labels has been generated.")

# Remove the temporary Froude and velocity columns
df_base.drop(columns=['tempF', 'tempV'], inplace=True)

# Plot 3: Boxplot of Manning's coefficient by slope percentage
plt.figure(figsize=(10, 8))
sns.boxplot(x='slope_perc', y='n.Manning', data=df_base, width=0.75, fliersize=1, palette='Set3')
sns.scatterplot(x='slope_perc', y='n.Manning', data=df_base, s=100, legend=False)

# Add text labels for flow rate
for i in range(df_base.shape[0]):
    plt.text(x=i % len(df_base['slope_perc'].unique()), y=df_base['n.Manning'][i], s=df_base['q_m3_h'][i],
             color='black', ha='right', va='bottom', fontsize=12, rotation=45)

plt.axhline(0, color='grey', linestyle='--', linewidth=0.75, alpha=0.6)  # Add a horizontal line at y=0
plt.title("Boxplot - Distribution of Manning's Coefficient", fontsize=16)
plt.xlabel("Slope (%)", fontsize=14)
plt.ylabel("Coefficient (-)", fontsize=14)
plt.xticks(rotation=0)
plt.grid(True)
plt.show()
print("Plot 3: Boxplot of Manning's Coefficient by Slope has been generated.")

# Plot 4: Boxplot of Kutter's coefficient by slope percentage
plt.figure(figsize=(10, 8))
sns.boxplot(x='slope_perc', y='m.Kutter', data=df_base, width=0.75, fliersize=1, palette='Set3')
sns.scatterplot(x='slope_perc', y='m.Kutter', data=df_base, s=100, legend=False)

plt.axhline(0, color='grey', linestyle='--', linewidth=0.75, alpha=0.6)
plt.title("Boxplot - Distribution of Kutter's Coefficient", fontsize=16)
plt.xlabel("Slope (%)", fontsize=14)
plt.ylabel("Coefficient (-)", fontsize=14)
plt.xticks(rotation=0)
plt.grid(True)
plt.show()
print("Plot 4: Boxplot of Kutter's Coefficient by Slope has been generated.")

# Plot 5: Boxplot of Bazin's coefficient by slope percentage
plt.figure(figsize=(10, 8))
sns.boxplot(x='slope_perc', y='s.Bazin', data=df_base, width=0.75, fliersize=1, palette='Set3')
sns.scatterplot(x='slope_perc', y='s.Bazin', data=df_base, s=100, legend=False)

plt.axhline(0, color='grey', linestyle='--', linewidth=0.75, alpha=0.6)
plt.title("Boxplot - Distribution of Bazin's Coefficient", fontsize=16)
plt.xlabel("Slope (%)", fontsize=14)
plt.ylabel("Coefficient (-)", fontsize=14)
plt.xticks(rotation=0)
plt.grid(True)
plt.show()
print("Plot 5: Boxplot of Bazin's Coefficient by Slope has been generated.")

# Create a new DataFrame for the normalized coefficients
df_norm = df_base[['slope_perc', 'n.Manning.norm', 'm.Kutter.norm', 's.Bazin.norm']]

# Convert data from "wide" format to "long" format for plotting
df_norm = df_norm.melt(id_vars=['slope_perc'], var_name='variable', value_name='value')

# Plot 6: Comparison of normalized coefficients (relative to the median)
plt.figure(figsize=(10, 8))
sns.boxplot(x='slope_perc', y='value', hue='variable', data=df_norm, palette='Set3', fliersize=3.5)

plt.axhline(0, color='grey', linestyle='--', linewidth=0.75, alpha=0.6)
plt.title("Boxplot - Comparison of Normalized Coefficients (Relative to Median)", fontsize=16)
plt.xlabel("Slope (%)", fontsize=14)
plt.ylabel("Normalized Coefficient (-)", fontsize=14)
plt.xticks(rotation=0)
plt.grid(True)
plt.show()
print("Plot 6: Comparison of Normalized Coefficients by Slope has been generated.")

# Round the relevant DataFrames to 3 decimal places
df_output = df_base.round(3)

print("\n/////////////////////////////////////////////////////////////")
print("END OF SCRIPT")
print("/////////////////////////////////////////////////////////////")
