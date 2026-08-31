import pandas as pd


# Load the CSV file into a DataFrame
csv_file = "Sample26_Absorbance_Spectrum.csv"
data = pd.read_csv(csv_file, low_memory=False)


# Drop any rows with NaN values
data.dropna(inplace=True)


# Initialize an empty 2D array
data_array = []


# Populate the 2D array with the DataFrame values
for i in range(min(201, len(data))):
   row = []
   for j in range(min(23, len(data.columns))):
       row.append(data.iloc[i, j])
   data_array.append(row)


# Print the 2D array
for row in data_array:
   print(row)


lowpeak = 0
highpeak = 0
max_val = 0.0


highpeak_sub = []
lowpeak_sub = []


# High Peak
for i in range(150, min(201, len(data_array))):
   for j in range(1, min(23, len(data_array[i]))):
       val = float(data_array[i][j])
       if val > max_val:
           max_val = val
           highpeak = i + 400  # 400 is added because of initial wavelength


# Ensure highpeak is within the valid range of indices
highpeak_index = highpeak - 400
if 0 <= highpeak_index < len(data_array):
   for i in range(1, min(23, len(data_array[highpeak_index]))):
       highpeak_cval = float(data_array[highpeak_index][i])
       lastcolumn_cval = float(data_array[200][i])
       highpeak_sub.append(highpeak_cval - lastcolumn_cval)


print(highpeak)
print(highpeak_sub)



