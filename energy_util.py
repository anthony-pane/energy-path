# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 09:13:46 2023

@author: Anthony J. Pane
"""

import numpy as np
import math
import scipy.interpolate as irp

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def continuous(energy_surface,method):
    """"Return interpolate object of energy surface"""
    unique_x = np.unique(energy_surface[:,0])
    unique_y = np.unique(energy_surface[:,1])
    energy_shape = (len(unique_x),
                    len(unique_y)
                    )
    energy = np.reshape(energy_surface[:,2],energy_shape)
    
    rgi = irp.RegularGridInterpolator((unique_x,unique_y),
                                      energy,method=method)
    return rgi
    

def value_energy(x_value,y_value,energy_surface,inter_func=None):

    x_unique = np.unique(energy_surface[:,0])
    y_unique = np.unique(energy_surface[:,1])
    x = find_nearest(x_unique, x_value)
    y = find_nearest(y_unique, y_value)
    
    index = np.where(( energy_surface[:,0] == x ) &
                     ( energy_surface[:,1] == y ) )[0][0]
    row = energy_surface[index]
    if inter_func is None:
        return row[2], row[3]
    else:
        return inter_func([x_value,y_value])[0], row[3]
    
    
def angstroms_to_degrees(angstrom,factor=0.9/.14):
    return angstrom * factor
    
def normalized_profile(profile,factor=0.9/.14,flip=False):
    """ normalize the profile to degree space """
    
    x = profile[0][0]
    y = angstroms_to_degrees(profile[0][1],factor)
    norm_profile = [[0.0,
                    profile[0][2],x,profile[0][1]]]
                    
    prev_x = x
    prev_y = y
    prev_distance = 0.0
    for position in profile[1:]:
        x = position[0]
        y = angstroms_to_degrees(position[1],factor)
        dx = (x - prev_x)**2
        dy = (y - prev_y)**2
        distance = np.sqrt(dx+dy) + prev_distance
        norm_profile.append([distance,
                             position[2],position[0],position[1]])
        prev_x = x
        prev_y = y
        prev_distance = distance
    
    norm_profile = np.asarray(norm_profile)
    if flip:
        norm_profile[:,0] = np.flip(norm_profile[:,0])
        
    return norm_profile

"""WARNING, CHECK IF FLIPPED"""
def step_coords(profile,norm_profile,factor,steps):
    total_distance = norm_profile[-1][0]
    step_distance = total_distance / ( steps - 2 )
    x = profile[0][0]
    y = profile[0][1]
    ener = profile[0][2]
    dist = 0
    coordinates = [[x,y,ener,dist]]
    for i,dist in enumerate(norm_profile[:,0]):    
        if dist >= step_distance:
            x = profile[i][0]
            y = profile[i][1]
            ener = profile[i][2]
            coordinates.append([x,y,ener,dist])
            step_distance += (total_distance / (steps - 1) )
    x = profile[-1][0]
    y = profile[-1][1]
    ener = profile[-1][2]
    coordinates.append([x,y,ener,dist])
    return np.asarray(coordinates)


