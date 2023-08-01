# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 14:01:02 2023

@author: Anthony J. Pane
"""
import numpy as np

class Spring:
    def __init__(self,k_constant,equilibrium_dist,dimension,power=2):
        
        self.type = "Spring"
        self.k = k_constant
        self.eq_dist = equilibrium_dist
        self.dim = dimension
        self.power = power
        
    def energy(self,node1,node2,power=2):
        if node1 is None:
            return 0.0
        if node2 is None:
            return 0.0
        if self.dim == 0:
            instant_dist = np.abs(node1.x - node2.x)
        elif self.dim == 1:
            instant_dist = np.abs(node1.y - node2.y)
            #print(instant_dist)
        return 0.5 * self.k * ((self.eq_dist - instant_dist)**self.power)
