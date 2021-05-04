#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--wt_file_name', type=str)
parser.add_argument('--tu_file_name', type=str)
parser.add_argument('--window_size', default = 10000, type=int)

args = parser.parse_args()

def read_depth_file(filename):
    f_wt_depth = open(filename, "r")

    chromosomes, positions, reads = [], [], []

    while True: 
        # Get next line from file
        line = f_wt_depth.readline()

        if not line:
            break
        (chrom, pos, read) = line.split("\t")
        chromosomes.append(chrom), positions.append(int(pos)), reads.append(int(read))

    f_wt_depth.close()
    return chromosomes, np.array(positions), np.array(reads)


def plot_read_depth(x, y, chrom, type_genome, window_size, color, im_name, ratio=False):
    plt.figure(figsize=(15, 10))
    plt.plot(x/1000000, y, 'o', markersize=5, markeredgewidth=0.7, markeredgecolor="w", color=color)

    plt.xlabel(f'Chromosome {chrom[0]} positions [Mb]')
    if not ratio:
        plt.ylabel('Read-depth') 
    else:
        plt.ylabel('Log2 ratio') 
    plt.title(f'Read depth {type_genome} (averaged every {window_size} bases)')
    plt.savefig(im_name)


def average_over_window(x, y, window_size):
    mod = len(x) % window_size

    x_mean = np.mean(x[:-mod].reshape(-1, window_size), axis=1)
    x_mean_rest = np.array(np.mean(x[-mod:]))[np.newaxis] 
    x_mean = np.concatenate((x_mean, x_mean_rest), axis=0)
    
    y_mean = np.mean(y[:-mod].reshape(-1, window_size), axis=1)
    y_mean_rest = np.array(np.mean(y[-mod:]))[np.newaxis] 
    y_mean = np.concatenate((y_mean, y_mean_rest), axis=0)
    return x_mean, y_mean


def main(args):
    wt_depth_filename = args.wt_file_name
    tu_depth_filename = args.tu_file_name
    window_size = args.window_size
    
    wt_chromosomes, wt_positions, wt_reads = read_depth_file(wt_depth_filename)
    tu_chromosomes, tu_positions, tu_reads = read_depth_file(tu_depth_filename)
    
    wt_mean_positions, wt_mean_reads = average_over_window(wt_positions, wt_reads, window_size)
    tu_mean_positions, tu_mean_reads = average_over_window(tu_positions, tu_reads, window_size)
    
    plot_read_depth(wt_mean_positions, wt_mean_reads, wt_chromosomes, "of normal genome", window_size, 'green', "normal_depth.png")
    plot_read_depth(tu_mean_positions, tu_mean_reads, tu_chromosomes, "of tumor genome", window_size, 'red', "tumor_depth.png")
    log_ratio = np.log2(tu_mean_reads/wt_mean_reads)
    plot_read_depth(tu_mean_positions, log_ratio, tu_chromosomes, "ratio", window_size, 'blue', "log2_ratio.png", ratio=True)


if __name__ == "__main__":
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)

