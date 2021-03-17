#!/usr/bin/env python
"""Final Project - Interactive Coordinate Match Module"""
"""Tzvetelina Dimitrova"""

"""Importing packages"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
from astropy import units as u
from astropy.coordinates import SkyCoord
import csv 
import pandas as pd
from astropy.coordinates import match_coordinates_sky
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

"""PART 1: Coordinate Matching Function"""
"""Function takes two datasets, an angle seperation, and outputs matched results in the first dataset"""
def match_stars(filename, filename2, sep):
    """Importing the datasets"""
    df = pd.read_csv(filename)
    RA = df["RA"]
    DEC = df["DEC"]
    data_rsgRADEC = pd.concat([RA,DEC], axis=1)

    df2 = pd.read_csv(filename2)
    RAstar = df2["RA"]
    DECstar = df2["DEC"]
    data2_rsgRADEC = pd.concat([RAstar,DECstar], axis=1)
    """Getting RA & DEC into skycoord"""
    rsg_coord = data_rsgRADEC.values.tolist()
    c1 = SkyCoord(rsg_coord, frame='icrs', unit= (u.hourangle, u.deg))
    
    rsg2_coord = data2_rsgRADEC.values.tolist()
    c2 = SkyCoord(rsg2_coord, frame='icrs', unit= (u.hourangle, u.deg))
    """Setting an angle seperation, finding matches"""
    max_sep = sep * u.arcsec
    idx, d2d, d3d = c2.match_to_catalog_sky(c1)
    sep_constraint = d2d < max_sep
    c2_matches = c2[sep_constraint]
    c1_matches = c1[idx[sep_constraint]]
    """Locating the matching stars in my dataset"""
    idx_matches = idx[sep_constraint]
    RSG_match = df.iloc[idx_matches]
    """Adding angle seperation and distance as columns"""
    angle_sep = d2d[sep_constraint]
    angle_sep_df = DataFrame(angle_sep, columns=['sep'])
    data_match = RSG_match.join(angle_sep_df.set_index(RSG_match.index), on=RSG_match.index)

    distance = d3d[sep_constraint]
    distance_df = DataFrame(distance, columns=['distance'])
    final_match = data_match.join(distance_df.set_index(RSG_match.index), on=RSG_match.index)
        
    return final_match
"""Make a new dataframe of the matched results"""
"""matched = match_stars('RSG1.csv','RSGpaper1.csv', 2.0)"""


"""PART 2: HR Diagram Function"""
"""Function takes original dataset, matches in the dataset, and data pointsizes, outputs HRD"""    
def hrd(filename, matches, pointsize,pointsize2):
    """Determining temperature & luminosity of my datasets stars"""
    data = pd.read_csv(filename)
    """Setting constant values"""
    R =  0.50
    reddening = 0.25
    distance_modulus = 5*np.log10(R * (10**6))-5
    K = data["K"]
    JK = data["J-K"]
    """Transforming colors to Bessel & Brett photometric system """
    JK_bb = (JK + 0.11) / 0.972
    K_bb = K + 0.044
    """De-reddening photometry"""
    JK_bb0 = (JK_bb) - 0.54* reddening
    K_bb0 = (K_bb) - 0.367*reddening
    M_K = K_bb0 - distance_modulus
    """Transforming to temperature, bolometric correction & magnitude"""
    temp = 5643.5 - (1807.1 * (JK_bb0))
    BC_K = 5.567 - (7.5686 * (10**(-4))) * temp
    M_bol = BC_K + M_K
    """Finding luminosity"""
    luminosity = (M_bol - 4.75) / -2.5
    T = np.log10(temp)
    L = luminosity
    
    """Determing temperature & luminosity of the matching stars"""    
    K_match = matches["K"]
    JK_match = matches["J-K"]
    JK_bb_match = (JK_match + 0.11) / 0.972
    K_bb_match = K_match + 0.044
    JK_bb0_match = (JK_bb_match) - 0.54* reddening
    K_bb0_match = (K_bb_match) - 0.367*reddening
    M_K_match = K_bb0_match - distance_modulus
    temp_match = 5643.5 - (1807.1 * (JK_bb0_match))
    BC_K_match = 5.567 - (7.5686 * (10**(-4))) * temp_match
    M_bol_match = BC_K_match + M_K_match
    luminosity_match = (M_bol_match - 4.75) / -2.5
    T_match = np.log10(temp_match)
    L_match = luminosity_match
    
    """Creating the HRD plot"""    
    fig, ax = plt.subplots(1,1)       
    fig.set_size_inches(13,8)
    ax.set_xlim(3.5,3.65)
    ax.set_ylim(3.9,5.5)
    ax.invert_xaxis()
    ax.set_xlabel("log(temperature)", fontsize=30)
    ax.set_ylabel("log(luminosity)", fontsize=30)
    plt.title('HRD of NGC6822 RSG Sample Matches', fontsize=30)

    ax.plot(T, L,color= "pink", marker="o",linestyle="None",markersize=pointsize, alpha=0.5);
    ax.plot(T_match,L_match,color="blue",marker="o",linestyle="None",markersize=pointsize2);

"""PART 3: Interactive HR Diagram"""
"""Function for interaction with dataset pointsizes"""
def interactive_hrd(filename, matches, pointsize, pointsize2):
    interact(hrd, filename = fixed(filename),matches = fixed(matches), pointsize=(0, 10, 1),pointsize2=(0, 10, 1));  
