import os
import shutil
import numpy as np
import scipy.io
import matplotlib.pyplot as plt


def reset_directory(dir_path, clean=True):
    """
    Ensure a directory exists at the specified path and optionally clean it.

    Args:
        dir_path (str): Path to the directory.
        clean (bool): If True, clear the directory if it exists;
        otherwise, ensure it exists without altering contents.
    """
    try:
        if clean and os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path, exist_ok=True)
    except Exception as e:
        print(f"Error resetting directory {dir_path}: {e}")


def read_timelog_our(filepath, num_networks):
    """
    Read time log entries from a file, structured by networks or blocks.

    Args:
        filepath (str): Path to the log file.
        num_networks (int): Number of networks (blocks) after each info line.

    Returns:
        list of dict: Each dictionary contains 'info' and
        'avg_computation_time' for a set of networks.
    """
    times = []
    try:
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            return []

        with open(filepath, 'r') as file:
            while True:
                line = file.readline().strip()
                if not line:
                    break
                time_block = [list(map(float, file.readline().strip().split()))
                              for _ in range(num_networks)]
                identifier_array = np.array(list(map(float, line.split())))
                time_array = np.array(time_block)
                entry = {
                    'info': identifier_array,
                    'avg_computation_time': np.mean(time_array[:, 1]),
                }
                times.append(entry)
    except Exception as e:
        print(f"Error reading or processing file {filepath}: {e}")
    return times


def read_timelog_cgbn(filepath):
    """
    Read average computation times from a MATLAB file.

    Args:
        filepath (str): Path to the MATLAB log file.

    Returns:
        numpy.ndarray: Array of average computation times.
    """
    try:
        times = scipy.io.loadmat(filepath)['average_times']
    except Exception as e:
        print(f"Error loading MATLAB file {filepath}: {e}")
        times = None
    return times


def read_timelog_bnbp(filepath):
    """
    Read time log from a file.

    Args:
        filepath (str): Path to the file.

    Returns:
        list of float: Times, convertedto float values.
    """
    times = []
    try:
        with open(filepath, 'r') as file:
            for line in file:
                line = line.strip()
                times = [float(value) for value in line.split(' ')]
                break
    except Exception as e:
        print(f"Error reading file {filepath} or processing data: {e}")
    return times


def plot_data(x_vals, y_vals, label, color, linestyle, marker):
    """
    Helper function to plot data with given styles.

    Args:
        x_vals (list): X-axis values.
        y_vals (list): Y-axis values.
        label (str): Label for the plot.
        color (str): Color of the plot line.
        linestyle (str): Line style of the plot.
        marker (str): Marker style of the plot.
    """
    plt.plot(x_vals, y_vals, label=label,
             color=color, linestyle=linestyle, marker=marker)


def visualize_combined(time_our, time_cgbn, time_bnbp, nodes, continuous,
                       target_dir='./Images/'):
    """
    Visualize computation times for various methods and configurations.

    Args:
        time_our (list of dict): Our method's time and info data.
        time_cgbn (numpy.ndarray): CGBN method's time data.
        time_bnbp (list): BayesNetBP method's computation times.
        nodes (list): List of nodes configurations.
        continuous (list): List of continuous variable proportions.
        target_dir (str): Directory to save the plots.
    """
    reset_directory(target_dir, clean=False)

    colors = ['#1f77b4', '#2ca02c', '#d62728']
    markers = ['o', '^', '*']
    line_styles = ['-', '--', ':']

    # Plot Our times
    for i, vertex in enumerate(nodes):
        x_vals, y_vals = [], []
        for entry in time_our:
            if entry['info'][0] == vertex and entry['info'][2] == 5:
                ratio = entry['info'][1]
                x_vals.append(ratio)
                y_vals.append(entry['avg_computation_time'])
        if x_vals:
            sorted_indices = np.argsort(x_vals)
            x_vals = np.array(x_vals)[sorted_indices]
            y_vals = np.array(y_vals)[sorted_indices]
            plot_data(x_vals, y_vals, 'Our Method',
                      colors[0], line_styles[0], markers[0])

    # Plot CGBN times
    for n in range(len(nodes)):
        times = [time_cgbn[n, p, 1] for p in range(len(continuous))]
        plot_data(continuous, times, 'CGBayesNets',
                  colors[1], line_styles[1], markers[1])

    # Plot BayesNetBP times
    plot_data(continuous, time_bnbp, 'BayesNetBP',
              colors[2], line_styles[2], markers[2])

    plt.xticks(continuous)
    plt.xlabel('Proportion of Continuous Variables')
    plt.ylabel('Computation Time (seconds, log scale)')
    plt.legend()
    plt.grid(False)
    plt.yscale('log')
    # plt.savefig(f'{target_dir}/comparison1.png', format='png', dpi=600)
    plt.savefig(f'{target_dir}/comparison1.eps', format='eps', dpi=600)
    plt.show()


if __name__ == "__main__":
    source1 = './BPMBN/BPMBN/results/'
    source2 = './CGBayesNets/'
    source3 = './BayesNetBP/'
    file_our = os.path.join(source1, 'timecost.txt')
    file_cgbn = os.path.join(source2, 'V50N1B100.mat')
    file_bnbp = os.path.join(source3, 'Time_V50E5.txt')
    num_networks = 100
    nodes = [50]
    continuous = [0, 0.25, 0.5, 0.75, 1]

    time_our = read_timelog_our(file_our, num_networks)
    time_cgbn = read_timelog_cgbn(file_cgbn)
    time_bnbp = read_timelog_bnbp(file_bnbp)

    visualize_combined(time_our, time_cgbn, time_bnbp, nodes, continuous)
