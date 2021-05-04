#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

#parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--wt_file_name', type=str)
parser.add_argument('--tu_file_name', type=str)
parser.add_argument('--window_size', default = 10000, type=int)
args = parser.parse_args()

# function to read the depth file line by line
# returns numpy arrays of chromosome types for each position, position index and the read-depth at the corresponding position
def read_depth_file(filename):
    f_depth = open(filename, "r")

    chromosomes, positions, reads = [], [], []
    while True: 
        # Get next line from file
        line = f_depth.readline()
        # break the while cycle at the end of the file
        if not line:
            break
        
        # split the line and append values to the lists
        (chrom, pos, read) = line.split("\t")
        chromosomes.append(chrom), positions.append(int(pos)), reads.append(int(read))

    f_depth.close()
    return chromosomes, np.array(positions), np.array(reads)

# function to plot the read depth averaged over a windows defined by window_size
# save the plot in PNG file with name defined in im_name
def plot_read_depth(x, y, chrom, type_genome, window_size, color, im_name, ratio=False):
    # plot settings - size of figure, labels, title
    plt.figure(figsize=(15, 10))
    plt.plot(x/1000000, y, 'o', markersize=5, markeredgewidth=0.7, markeredgecolor="w", color=color)

    plt.xlabel(f'Chromosome {chrom[0]} positions [Mb]')
    if not ratio:
        plt.ylabel('Read-depth') 
    else:
        plt.ylabel('Log2 ratio') 
    plt.title(f'Read depth {type_genome} (averaged every {window_size} bases)')
    plt.savefig(im_name)

# function to take average of the reads (and positions) over the window - shoul be rigid against window size, which does not divide the length of the data
def average_over_window(x, y, window_size):
    # take modulo to deal with window size, which do not divide the length of the data
    mod = len(x) % window_size
    
    # take the reads without the values in the last window (if neccesary), get means of the windows
    # concatenate with mean of the last window
    if mod == 0:
        x_mean = np.mean(x.reshape(-1, window_size), axis=1)   
    else:
        x_mean = np.mean(x[:-mod].reshape(-1, window_size), axis=1)
        x_mean_rest = np.array(np.mean(x[-mod:]))[np.newaxis] 
        x_mean = np.concatenate((x_mean, x_mean_rest), axis=0)
    
    if mod == 0:
        y_mean = np.mean(y.reshape(-1, window_size), axis=1)
    else:
        y_mean = np.mean(y[:-mod].reshape(-1, window_size), axis=1)
        y_mean_rest = np.array(np.mean(y[-mod:]))[np.newaxis] 
        y_mean = np.concatenate((y_mean, y_mean_rest), axis=0)
    return x_mean, y_mean

# call all methods to produce the results
def main(args):
    # assign arguments
    wt_depth_filename = args.wt_file_name
    tu_depth_filename = args.tu_file_name
    window_size = args.window_size
    
    # read files into three numpy arrays
    wt_chromosomes, wt_positions, wt_reads = read_depth_file(wt_depth_filename)
    tu_chromosomes, tu_positions, tu_reads = read_depth_file(tu_depth_filename)
    
    # average the reads and positions over the windows 
    wt_mean_positions, wt_mean_reads = average_over_window(wt_positions, wt_reads, window_size)
    tu_mean_positions, tu_mean_reads = average_over_window(tu_positions, tu_reads, window_size)
    
    # calculate log2 ratio, not dividing by zero
    log_ratio = np.log2(np.divide(tu_mean_reads, wt_mean_reads, where=wt_mean_reads!=0))
    
    # plot the data and save the final plots
    plot_read_depth(wt_mean_positions, wt_mean_reads, wt_chromosomes, "of normal genome", window_size, 'green', "normal_depth.png")
    plot_read_depth(tu_mean_positions, tu_mean_reads, tu_chromosomes, "of tumor genome", window_size, 'red', "tumor_depth.png")
    plot_read_depth(tu_mean_positions, log_ratio, tu_chromosomes, "ratio", window_size, 'blue', "log2_ratio.png", ratio=True)


if __name__ == "__main__":
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)

