import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_csv(file_path):
   """Load the CSV file into a DataFrame."""
   return pd.read_csv(file_path, low_memory=False)


def get_concentration_info():
   """Get the original concentration and stock solution percentage from the user."""
   ogCon = float(input("What was the original concentration of the substance in moles:  "))
   ogPer = float(input("Enter percentage of that concentration in your stock solution? Do not include the percent symbol: "))
   per = ogPer / 100.0
   conInMol = per * ogCon
   conInMmol = conInMol * 1000
   print(conInMmol)
   return conInMmol


def preprocess_data(data):
   """Drop rows with NaN values and convert the DataFrame to a 2D array."""
   data.dropna(inplace=True)
   data_array = [
       [data.iloc[i, j] for j in range(min(23, len(data.columns)))]
       for i in range(min(201, len(data)))
   ]
   return data_array


def find_high_peak(data_array):
   """Find the high peak value and its corresponding row."""
   max_val = 0.0
   highpeak = 0
   for i in range(1, min(201, len(data_array))):
       val = float(data_array[i][1])
       if val > max_val:
           max_val = val
           highpeak = i + 400  # 400 is added because of initial wavelength
   print("The high peak is: ", highpeak)
   return highpeak


def calculate_peak_subtraction(data_array, peak, last_column_index=200):
   """Calculate the peak subtraction values."""
   peak_sub = []
   peak_index = peak - 400
   if 0 <= peak_index < len(data_array):
       for i in range(1, min(23, len(data_array[peak_index]))):
           peak_cval = float(data_array[peak_index][i])
           lastcolumn_cval = float(data_array[last_column_index][i])
           peak_sub.append(peak_cval - lastcolumn_cval)
   return peak_sub


def find_low_peak(data_array):
   """Find the low peak value and its corresponding row."""
   max_val = 0.0
   lowpeak = 0
   for i in range(2, min(131, len(data_array))):
       for j in range(1, min(23, len(data_array[i]))):
           val = float(data_array[i][j])
           if val > max_val:
               max_val = val
               lowpeak = i + 400  # 400 is added because initial wavelength
   print("The low peak is: ", lowpeak)
   return lowpeak


def calculate_absorbance(highpeak_sub, lowpeak_sub):
   """Calculate the absorbance values."""
   absorbance = [highpeak_sub[i] / lowpeak_sub[i] for i in range(len(lowpeak_sub))]
   return absorbance


def perform_linear_regression(absorbance):
   """Perform linear regression on the last three absorbance values."""
   x = np.array([0, 5, 10])  # Assuming x values for the last three measurements
   y = np.array(absorbance[-3:])  # Last three absorbance values
   coefficients = np.polyfit(x, y, 1)  # Linear fit (degree 1)
   slope, intercept = coefficients
   print(f"Linear Regression: slope = {slope}, intercept = {intercept}")
   return slope, intercept


def main(csv_file):
   data = load_csv(csv_file)
   conInMmol = get_concentration_info()
   data_array = preprocess_data(data)
   highpeak = find_high_peak(data_array)
   highpeak_sub = calculate_peak_subtraction(data_array, highpeak)
   #print(highpeak_sub)
   #print("\n")
  
   lowpeak = find_low_peak(data_array)
   lowpeak_sub = calculate_peak_subtraction(data_array, lowpeak)
   #print(lowpeak_sub)
   #print("\n")
  
   absorbance = calculate_absorbance(highpeak_sub, lowpeak_sub)
   print(absorbance)
  
   # Perform linear regression on the last three absorbance values
   perform_linear_regression(absorbance)


# Usage
csv_file = input("Enter file name: ")
main(csv_file)
