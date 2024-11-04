import sys
import numpy as np
import csv

def calculate_statistics(filename):
    times = []

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        for row in reader:
            try:
                print(row)
                times.append(float(row[0].strip()))
            except ValueError:
                continue
    print(times)
    # Ensure there are enough values to drop the highest and lowest
    if len(times) < 3:
        raise ValueError("Not enough data points to calculate statistics after dropping extremes.")

    # Drop the highest and lowest times
    times.sort()
    trimmed_times = times[1:-1]  # Remove the first (lowest) and last (highest)

    return np.mean(trimmed_times), np.std(trimmed_times)

def main(log_prefix):
    log_filename = f"{log_prefix}.csv"
    try:
        mean_time, std_dev_time = calculate_statistics(log_filename)
        print(f"Mean time: {mean_time:.6f}")
        print(f"Standard Deviation: {std_dev_time:.6f}")
    except FileNotFoundError:
        print(f"Error: The file {log_filename} does not exist.")
    except ValueError as ve:
        print(f"Error: {ve}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_stats.py <log_prefix>")
    else:
        main(sys.argv[1])
