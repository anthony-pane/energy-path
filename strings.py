# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 17:36:00 2023

@author: Anthony J. Pane
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import nodes
import springs
from energy_util import *
from scipy.interpolate import CubicSpline

class String:
    def __init__(self,energy_surface,positions,inter_function=None):
        """
        inter_function is either “linear”, 
            “nearest”, “slinear”, “cubic”, 
            “quintic” or “pchip”
        """
        self.energy_surface = energy_surface
        self.positions = positions
        if inter_function is not None:
            self.inter_function = continuous(self.energy_surface,
                                             inter_function)
        else:
            self.inter_function = None
        self.node_list = nodes.Node(self.energy_surface,
                                 self.positions,prev_node=None,
                                 node_index=0, inter_function=None)
        self.energies = []
        self.dx = 0
        self.dy = 0

    
    def find_node(self,index):
        """ Return the node at an index in node list """
        node = self.node_list
        try:
            while(node is not None):
                if node.index == index:
                    return node
                node = node.next_node
        except AttributeError:
            print("Index not in node list")
            
    
    def node_pmf_energy(self,index):
        """ Return the pmf energy of node at index """
        return self.find_node(index).pmf_energy
    
    def node_neighbors(self,index):
        """ Return all of the neighbors of the node at index """
        return self.find_node(index).neighborhood(self.dx,
                                                  self.dy,
                                                  periodic_x=False,
                                                  periodic_y=False)
    
    def node_spring_energy(self,index):
        """ Return the energy of the spring(s) at node index """
        node = self.find_node(index)
        spring_energy = 0
        for spring in self.energies:
            if index > 0:
                prev_node = self.find_node(index-1)
                spring_energy += spring.energy(prev_node,node)
            if index < len(self.positions) - 1:
                next_node = self.find_node(index+1)
                spring_energy += spring.energy(node,next_node)  
        return spring_energy

    def find_minimum(self,node):
        """
        Minimize the node to a local minimum. Local is x+dx
        and y+dy. The energy of the pmf and any additional springs
        are included.
        """
        original_x = node.x
        original_y = node.y
        min_x = node.x
        min_y = node.y
        min_energy = node.pmf_energy + node.external_energy
        for neighbor in node.neighborhood(self.dx,self.dy):
            node.x = neighbor[0]
            node.y = neighbor[1]
            node.update_energy(self.energy_surface,self.inter_function)
            node_energy = node.pmf_energy + node.external_energy
            if node_energy < min_energy:
                min_x = node.x
                min_y = node.y
                min_energy = node_energy
            node.x = original_x
            node.y = original_y
            node.update_energy(self.energy_surface,self.inter_function)
        return min_x,min_y
    
    def move_node(self,node,newx,newy):
        if not node.static:
            node.x = newx
            node.y = newy
            node.update_energy(self.energy_surface,self.inter_function)
    
    def minimize(self,nsteps):
        """
        Move each node to a minimum energy according to the PMF
        and the added spring parameters. Each node finds the lowest
        energy in x+dy and y+dy. Repeat for nsteps.
        """
        for step in range(nsteps):
            node = self.node_list
            new_positions = []
            while(node is not None):
                new_x,new_y = self.find_minimum(node)
                new_positions.append([new_x,new_y])
                node = node.next_node
                if node is None:
                    break
            node = self.node_list
            while(node is not None):
                self.move_node(node,
                               new_positions[node.index][0],
                               new_positions[node.index][1])
                node = node.next_node
                if node is None:
                    break
        self.update_positions()
    
    def pmf_energy(self):
        """ 
        Return the total energy contribution from the energy surface
        """
        total_energy = 0
        node = self.node_list
        while(node is not None):
            total_energy += node.pmf_energy
            node = node.next_node
        return total_energy

    def external_energy(self):
        """
        Return the external energy contribution 
        on all of the nodes.
        """
        total_external = 0
        node = self.node_list
        while(node is not None):
            total_external += node.external_energy
            node = node.next_node
            if node is None:
                break
        return total_external
    
    
    def add_energy(self,external_energy):
        self.energies.append(external_energy)
        
        
    def update_external(self):
        node = self.node_list
        node.external_energies = []
        for ener in self.energies:
            node.add_external_energies(ener)
        node.update_energy(self.energy_surface,self.inter_function)
        while(node is not None):
            node = node.next_node
            if node is None:
                break
            node.external_energies = []
            for ener in self.energies:
                node.add_external_energies(ener)
            node.update_energy(self.energy_surface,self.inter_function)
    
    
    def update_positions(self):
        new_positions = []
        node = self.node_list
        while(node is not None):
            new_positions.append([node.x,node.y])
            node = node.next_node
            if node is None:
                break
        self.positions = np.asarray(new_positions)
    
        
    def node_profile(self,nodes=-1):
        """Return the x,y,energy,probability for each node"""
        profile = []
        node =self.node_list
        while(node is not None):
            profile.append([node.x,node.y,
                            node.pmf_energy,node.external_energy,
                            node.probability])
            node = node.next_node
            if node is None:
                break
        profile = np.asarray(profile)
        return profile

            
            
    """
    Also include methods to remove external forces 
    """
    
    def set_dx(self,dx):
        self.dx = dx

    def set_dy(self,dy):
        self.dy = dy
        
    def plot(self,title=' '):
        """
        Plot a contour of the energy surface with the band.
        """      
        origin = 'lower'
        minlevel = 0
        maxlevel = 12
        steplevel = 1.0
        levels = np.arange(minlevel, maxlevel + steplevel, steplevel)   
        fig1, ax2 = plt.subplots(layout='constrained')
        CS = ax2.tricontourf(self.energy_surface[:,0],
                       self.energy_surface[:,1],
                       self.energy_surface[:,2],
                       vmin=0,
                       vmax=12,
                       levels = levels,
                       origin=origin)
        CS2 = ax2.tricontour(CS, levels=levels, 
                             colors = 'black',
                             origin=origin)
        ax2.set_title(title)
        cbar = fig1.colorbar(CS)
        cbar.add_lines(CS2)
        #print(np.shape(self.positions))
        plt.plot(self.positions[:,0],self.positions[:,1],'bo')
            
            
            
            