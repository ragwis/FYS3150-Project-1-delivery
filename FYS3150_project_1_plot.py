# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 19:52:44 2020

@author: rawis
"""

#import numpy as np
import matplotlib.pyplot as plt
#import sys


def extract_data(filename):
    """stores data from file in arrays for u, v and x"""
    infile = open(filename, 'r') #open and read file
    n = int(infile.readline())
    log10_h = float(infile.readline())
    x = [float(i) for i in infile.readline().split()]
    u_exact =   [float(i) for i in infile.readline().split()]
    u_standard =  [float(i) for i in infile.readline().split()]
    relative_error_standard= float(infile.readline())
    u_specific =   [float(i) for i in infile.readline().split()]
    relative_error_specific= float(infile.readline())
    #u_LU =   [float(i) for i in infile.readline().split()]
    #relative_error_LU= float(infile.readline())
    infile.close()
    #return x, u_exact, u_standard, u_specific, u_LU, relative_error_standard, relative_error_specific, relative_error_LU, n, log10_h
    return x, u_exact, u_standard, u_specific, relative_error_standard, relative_error_specific, n, log10_h

#filename = sys.argv[1]
def plot_data(filename):
    x, u_exact, u_standard, u_specific, relative_error_standard, relative_error_specific, n, log10_h = extract_data(filename)
    plt.plot(x, u_standard, label='$v_{general}$')
    #print(u_specific)
    plt.plot(x, u_specific, label='$v_{specific}$')
    #plt.plot(x, u_LU, label='$u_{computed_{LU}}$')
    plt.plot(x, u_exact, label='$u_{exact}$')
    plt.title(f'n = {n}')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.legend()
    plt.show()
    
def plot_error(filename):
    x, u_exact, u_standard, u_specific, relative_error_standard, relative_error_specific, n, log10_h = extract_data(filename)
    plt.scatter(log10_h, relative_error_standard)
    plt.xlabel('$\log_{10}(h)$')
    plt.ylabel('$\epsilon$')
    plt.show()
    print(f'log10_h = {log10_h}')
    print(f'relative_error_standard = {relative_error_standard}')

ns = [10**i for i in range(1,8)]
filenames = [ "data_n="+str(n)+".txt" for n in ns]
#for filename in filenames:
    #plot_error(filename)

plot_data(filenames[2])

