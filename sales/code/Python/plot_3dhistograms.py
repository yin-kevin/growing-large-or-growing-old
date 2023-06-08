# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:48:01 2022

@author: ky297
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# directory
directory_path = "D:\Kevin\Firm_Growth_Project\_Final\Sales\Output"
os.chdir(directory_path)

# import growth rate data
df_growth_rates = pd.read_csv(r"D:\Kevin\Firm_Growth_Project\_Final\Sales\Data\clean\new_vat_data_cleaned.csv")
df_growth_rates_cont = df_growth_rates[['year', 'gro_ihs',  'gro_halt',  'gro_simpl', 'gro_log', 'entry', 'exit']]
df_growth_rates = df_growth_rates[['year', 'gro_ihs',  'gro_halt',  'gro_simpl', 'gro_log']]


#%% all firms (including entrants and exiting)
for r in ['ihs', 'halt']:
    
    # determine bounds based on the type of growth rate
    if r == 'ihs':
        low = -20
    elif r == 'halt':
        low = -2
    elif r == 'simpl':
        low = -1
    elif r == 'log':
        low = -5
    
    
    
    # make the bounds symmetric across 0
    high = -low
    
    # construct histogram data
    hist, xedges, yedges = np.histogram2d(df_growth_rates['gro_{}'.format(r)], 
                                          df_growth_rates['year'], 
                                          bins=[50, 22], 
                                          range=[[low, high],[1996, 2018]])
    
    # Construct arrays for the anchor positions of the bars.
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    dz = hist.ravel()
    zpos = np.zeros(len(dz))
     
    # bar width and length
    dx = xedges[1]-xedges[0]
    dy = yedges[1]-yedges[0]
    
    # Construct arrays with the dimensions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='lightskyblue', edgecolor='black', lw=0.3)
    ax.view_init(elev=15, azim=-70)
    
    # ticksize
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    
    # location of tick labels
    ax.tick_params(axis="y",direction="in", pad=-3)
    ax.tick_params(axis="y",direction="in", pad=-3)
    ax.tick_params(axis="x",direction="in", pad=-3)
    ax.tick_params(axis='z',labelsize=6)
    
    # turn the y-ticks into integers from floats and put them closer
    y_ticks = [int(i) for i in ax.get_yticks()]
    
        
    ax.set_yticklabels(y_ticks, 
                       verticalalignment='baseline',
                       horizontalalignment='left')
    
    # title and axis labels
    ax.set_title("frequency of {} growth-rates".format(r), fontsize=10, y=0.97, x=0.51)
    ax.set_xlabel("growth rate", fontsize=7)
    ax.xaxis.set_rotate_label(False)
    ax.xaxis.labelpad=-4
    ax.set_ylabel("year", fontsize=7)
    ax.yaxis.labelpad=0
    
    plt.show()
    fig.savefig(r"D:\Kevin\Firm_Growth_Project\_Final\Employment\Output\figures\Histogram\histogram_3D_{}.png".format(r), bbox_inches='tight', dpi=300)
    


#%% no entrants, no exiting firms
df_growth_rates_cont = df_growth_rates_cont[(df_growth_rates_cont['entry'] == 0) & (df_growth_rates_cont['exit'] == 0)]

for r in ['ihs', 'halt', 'simpl', 'log']:
    
    # determine bounds based on the type of growth rate
    if r == 'ihs':
        low = -20
    elif r == 'halt':
        low = -2
    elif r == 'simpl':
        low = -1
    elif r == 'log':
        low = -5
    
    
    
    # make the bounds symmetric across 0
    high = -low
    
    # construct histogram data
    hist, xedges, yedges = np.histogram2d(df_growth_rates_cont['gro_{}'.format(r)], 
                                          df_growth_rates_cont['year'], 
                                          bins=[50, 22], 
                                          range=[[low, high],[1996, 2018]])
    
    # Construct arrays for the anchor positions of the bars.
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    dz = hist.ravel()
    zpos = np.zeros(len(dz))
     
    # bar width and length
    dx = xedges[1]-xedges[0]
    dy = yedges[1]-yedges[0]
    
    # Construct arrays with the dimensions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='lightskyblue', edgecolor='black', lw=0.3)
    ax.view_init(elev=15, azim=-70)
    
    # ticksize
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    
    # location of tick labels
    ax.tick_params(axis="y",direction="in", pad=-3)
    ax.tick_params(axis="y",direction="in", pad=-3)
    ax.tick_params(axis="x",direction="in", pad=-3)
    ax.tick_params(axis='z',labelsize=6)
    
    # turn the y-ticks into integers from floats and put them closer
    y_ticks = [int(i) for i in ax.get_yticks()]
    
        
    ax.set_yticklabels(y_ticks, 
                       verticalalignment='baseline',
                       horizontalalignment='left')
    
    # title and axis labels
    ax.set_title("frequency of {} growth-rates".format(r), fontsize=10, y=0.97, x=0.51)
    ax.set_xlabel("growth rate", fontsize=7)
    ax.xaxis.set_rotate_label(False)
    ax.xaxis.labelpad=-4
    ax.set_ylabel("year", fontsize=7)
    ax.yaxis.labelpad=0
    
    plt.show()
    fig.savefig(r"D:\Kevin\Firm_Growth_Project\_Final\Employment\Output\figures\Histogram\histogram_3D_{}_cont.png".format(r), bbox_inches='tight', dpi=300)
    

