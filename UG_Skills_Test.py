import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

file_Path = "data.csv"


# Load dataset from data.csv




# #1 Load file
def load_Dataset(file_Path):
    with open(file_Path, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  
        data = list(reader)
    return header, data

header, data = load_Dataset('data.csv')

# #2 Split dataset
def split_dataset(data):
    right_Eye_x = [float(row[1]) for row in data if row[0].startswith('Right_Calibration')]
    right_Eye_y = [float(row[2]) for row in data if row[0].startswith('Right_Calibration')]
    left_Eye_x = [float(row[3]) for row in data if row[0].startswith('Left_Calibration')]
    left_Eye_y = [float(row[4]) for row in data if row[0].startswith('Left_Calibration')]
    time = [float(row[5]) for row in data if row[0].startswith('Test_Data')]
    return right_Eye_x, right_Eye_y, left_Eye_x, left_Eye_y, time

# Split dataset
right_Eye_x, right_Eye_y, left_Eye_x, left_Eye_y, time = split_dataset(data)

# #3 Generate calibration coefficients
def generate_Calibration_Coefficients(x_Calibrations):
    return np.mean(x_Calibrations)

calibration_Coefficient_right = generate_Calibration_Coefficients(right_Eye_x)
calibration_Coefficient_left = generate_Calibration_Coefficients(left_Eye_x)

# #4 Butterworth filter
def butterworth_filter(data, order, cutoff_Frequency, sampling_Frequency):
    nyquist = 0.5 * sampling_Frequency
    normal_Cutoff = cutoff_Frequency / nyquist
    b, a = butter(order, normal_Cutoff, btype='low', analog=False)
    return filtfilt(b, a, data, padlen=2)
  # Filter each eye's data using a Butterworth filter
order = 3
cutoff_Frequency = 1
sampling_Frequency = 70

right_Eye_Filtered = butterworth_filter(right_Eye_x, order, cutoff_Frequency, sampling_Frequency)
left_Eye_Filtered = butterworth_filter(left_Eye_x, order, cutoff_Frequency, sampling_Frequency)

# #5 Center data around 0
def center_data(data):
    return data - np.mean(data)
  
right_Eye_centered = center_data(right_Eye_x)
left_Eye_centered = center_data(left_Eye_x)

# #6 Calibrate data using coefficients
def calibrate_Data(data, calibration_Coefficient):
    return data / calibration_Coefficient
 
right_Eye_calibrated = calibrate_Data(right_Eye_x, calibration_Coefficient_right)
left_Eye_calibrated = calibrate_Data(left_Eye_x, calibration_Coefficient_left)

# #7 Calculate final data
def calculate_final_Data(right_Eye_data, left_Eye_data):
    return right_Eye_data - left_Eye_data
  # Calculate final data using right eye - left eye
final_Data = calculate_final_Data(right_Eye_calibrated, left_Eye_calibrated)

# #8 Generate plots
def generate_plots(right_Eye_raw, left_Eye_raw, right_Eye_filtered, left_Eye_filtered,
                   right_Eye_centered, left_Eye_centered, right_Eye_calibrated, left_Eye_calibrated,
                   final_Data):
    # Figure 1
    plt.figure(1, figsize=(10, 8))

    # Subplot 1: Raw data for each eye
    plt.subplot(4, 1, 1)
    plt.plot(right_Eye_raw, color='red', label='Right Eye')
    plt.plot(left_Eye_raw, color='blue', label='Left Eye')
    plt.title('Raw Data')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.legend()

    # Subplot 2: Filtered data for each eye
    plt.subplot(4, 1, 2)
    plt.plot(right_Eye_filtered, color='red', label='Right Eye')
    plt.plot(left_Eye_filtered, color='blue', label='Left Eye')
    plt.title('Filtered Data')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.legend()

    # Subplot 3: Centered data for each eye
    plt.subplot(4, 1, 3)
    plt.plot(right_Eye_centered, color='red', label='Right Eye')
    plt.plot(left_Eye_centered, color='blue', label='Left Eye')
    plt.title('Centered Data')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.legend()

    # Subplot 4: Calibrated data for each eye
    plt.subplot(4, 1, 4)
    plt.plot(right_Eye_calibrated, color='red', label='Right Eye')
    plt.plot(left_Eye_calibrated, color='blue', label='Left Eye')
    plt.title('Calibrated Data')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.legend()

    plt.tight_layout()

    # Figure 2
    plt.figure(2, figsize=(10, 8))
    plt.plot(final_Data, label='Final Data', color='black')
    plt.title('Final Data')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.legend()

    plt.show()

# Generate plots
generate_plots(right_Eye_x, left_Eye_x, right_Eye_Filtered, left_Eye_Filtered, right_Eye_centered, left_Eye_centered, right_Eye_calibrated, left_Eye_calibrated,final_Data)
