# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 10:54:16 2023

@author: Anthony J. Pane
"""
import numpy as np
from dijkstra import *
from energy_util import *
import scipy.interpolate as irp


class Node:
    def __init__(self,energy_surface,positions,
                 prev_node,node_index,
                 inter_function=None):
        
        # static first and last node
        if node_index == len(positions) - 1:
            self.next_node = None
            self.prev_node = prev_node
            self.static = True
        elif node_index == 0:
            self.next_node = Node(energy_surface,positions,self,node_index+1)
            self.prev_node = None
            self.static = True
        else:
            self.next_node = Node(energy_surface,positions,self,node_index+1)
            self.prev_node = prev_node
            self.static = False
        
        #self.energy_surface = energy_surface
        self.x_min = np.min(energy_surface[:,0])
        self.x_max = np.max(energy_surface[:,0])
        self.y_min = np.min(energy_surface[:,1])
        self.y_max = np.max(energy_surface[:,1])
        
        
        self.x = positions[node_index][0]
        self.y = positions[node_index][1]
        self.index = node_index
        self.pmf_energy,self.probability = value_energy(self.x,
                                                        self.y,
                                                        energy_surface,
                                                        inter_function)
        self.external_energies = []
        self.external_energy = 0
        
        
        
    def update_energy(self,energy_surface,inter_function):
        self.pmf_energy,self.probability = value_energy(self.x,
                                                        self.y,
                                                        energy_surface,
                                                        inter_function)
        self.external_energy = 0
        for spring in self.external_energies:
            self.external_energy += spring.energy(self.prev_node,self)
            self.external_energy += spring.energy(self,self.next_node)
        
    def add_external_energies(self,external_energy):
        self.external_energies.append(external_energy)
        #print(self.external_energies[0].energy(self,self.next_node))
        """
        spring_energy = 0
        for spring in self.external_energies:
            if index > 0:
                spring_energy += spring.energy(self.prev_node,self)
            if index < len(self.positions) - 1:
                spring_energy += spring.energy(node,self.next_node)  
        return spring_energy
        """
        
    def neighborhood(self,dx,dy,periodic_x=False,periodic_y=False):
        """
        Return all of the
        neighboring x and y indexes 
        as an array of [x_index,y_index].
        Fix periodicity at a later date.
        """
        neighbors = []
        try:
            xrange = [x_val + self.x for x_val in [-1*dx,0,dx]]
            yrange = [y_val + self.y for y_val in [-1*dy,0,dy]]
        except:
            print(self.x)
            print(self.y)
        
        for xval in xrange:
            if xval >= self.x_min and xval <= self.x_max: # change here for periodic-x
                for yval in yrange:
                    if yval >= self.y_min and yval <= self.y_max: # here for periodic-y
                        neighbors.append([xval,yval])
        return neighbors                        
