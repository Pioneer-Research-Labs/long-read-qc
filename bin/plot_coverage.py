#!/usr/bin/env python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import subprocess



def plot_depth(depth_report, output_name, title, genome_size):
    """Plot genome Depth across genome.

    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        title (str): Title of the plot.
        genome_size (int): Genome size.

    """
    data = load_depth(depth_report, genome_size)
    y_label = "Count of Inserts"
    plt.subplots(figsize=(12, 6))

    plt.title(title)


    sns_plot = sns.lineplot(x=range(len(data)), y=data)
    sns_plot.set(xlabel='Genome Position (bp)', ylabel=y_label)

    plt.savefig(output_name, bbox_inches='tight')
    plt.close()


def load_depth(depth_report, rows):
    """Load depth report.

    Args:
        depth_report (str): Path to samtool's depth file.

    Returns:
        list: List with depth.

    """
    no_rows = rows

    a = np.zeros((no_rows), dtype=np.int64)

    with open(depth_report) as f:
        for i, line in enumerate(f):
            pos, count = line.rstrip().split("\t")[1:]
            a[i] =  int(count)
    return a.tolist()

if __name__ == '__main__':

    command_args = ['wc', '-l', sys.argv[1]]

    process = subprocess.Popen(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, errors = process.communicate()
    count = int(output.split()[0])
    plot_depth(sys.argv[1], sys.argv[2], sys.argv[3], count))
    plt.close()
