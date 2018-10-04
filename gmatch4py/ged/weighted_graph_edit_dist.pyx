# -*- coding: UTF-8 -*-

import sys

import networkx as nx
import numpy as np
cimport numpy as np
from .abstract_graph_edit_dist cimport AbstractGraphEditDistance
from ..base cimport intersection,union_



cdef class WeightedGraphEditDistance(AbstractGraphEditDistance):

    def __init__(self,node_del,node_ins,edge_del,edge_ins, edge_weight_coef=1, weighted=False):
        AbstractGraphEditDistance.__init__(self,node_del,node_ins,edge_del,edge_ins)
        self.edge_weight_coef=edge_weight_coef
        self.weighted=weighted

    cpdef double substitute_cost(self, node1, node2, G, H):
        return self.relabel_cost(node1, node2, G, H)

    def add_edges(self,node1,node2,G):
        R=nx.create_empty_copy(G)
        try:
            R.add_edges_from(G.edges([node1,node2], data=True))
        except Exception as e:
            # To counter bug with a None for attribute... weird ??
            arr_=G.edges([node1,node2], data=True) + G.in_edges([node1,node2], data=True)
            new_list=[]
            for item in arr_:
                new_list.append((item[0],item[1], item[2]))
            R.add_edges_from(new_list)
        return R

    def relabel_cost(self, node1, node2, G, H):
        # Check if nodes match both in terms of id and degree
        if (node1 == node2 
            and G.degree(node1, weight=("weight" if self.weighted else None))
                    == H.degree(node2, weight=("weight" if self.weighted else None))):
            # we need to also check if there is an edit cost for the matching edges
            g_out_weight = sum([edge[-1]['weight'] for edge in G.edges(node1, data=True)])
            h_out_weight = sum([edge[-1]['weight'] for edge in H.edges(node2, data=True)])
            weight_out_diff = abs(g_out_weight - h_out_weight)
            return weight_out_diff*self.edge_weight_coef
        # If node ids are the same but degrees differ, we need to take into account the missing edges' costs        
        elif (node1 == node2 
            and G.degree(node1, weight=("weight" if self.weighted else None)) 
                != H.degree(node2, weight=("weight" if self.weighted else None))):
            R = self.add_edges(node1,node2,G)
            R2 = self.add_edges(node1,node2,H)
            inter_=intersection(R,R2).number_of_edges()
            add_diff=abs(R2.number_of_edges()-inter_)
            del_diff=abs(R.number_of_edges()-inter_)
            weight_r = sum([edge[-1]['weight'] for edge in R.edges(data=True)])
            weight_r2 = sum([edge[-1]['weight'] for edge in R2.edges(data=True)])
            weight_diff = abs(weight_r - weight_r2)
            return (add_diff*self.edge_ins)+(del_diff*self.edge_del)+(weight_diff*self.edge_weight_coef)


        # For each missing edge we penalize both for substitutions and missing cost. 
        if (node1,node2) in G.edges():
            return self.node_ins+self.node_del \
                        +G.get_edge_data(node1, node2)['weight']*self.edge_weight_coef
        if (node2,node1) in G.edges():
            return self.node_ins+self.node_del \
                        +G.get_edge_data(node2, node1)['weight']*self.edge_weight_coef
        if not node2 in G:
            nodesH=list(H.nodes())
            index=nodesH.index(node2)
            return self.node_del+self.node_ins+self.insert_cost(index,index,nodesH,H)
        return sys.maxsize

    cdef double delete_cost(self, int i, int j, nodesG, G):
        if i == j:
            edge_weight = sum([edge[-1]['weight'] for edge in G.edges(nodesG[i], data=True)]) \
                            + sum([edge[-1]['weight'] for edge in G.in_edges(nodesG[i], data=True)])
            return self.node_del+(G.degree(nodesG[i],weight=("weight" if self.weighted else None))*self.edge_del)+edge_weight*self.edge_weight_coef
        return sys.maxsize

    cdef double insert_cost(self, int i, int j, nodesH, H):
        if i == j:
            edge_weight = sum([edge[-1]['weight'] for edge in H.edges(nodesH[j], data=True)]) \
                            + sum([edge[-1]['weight'] for edge in H.in_edges(nodesH[j], data=True)])
            deg=H.degree(nodesH[j],weight=("weight" if self.weighted else None))
            if isinstance(deg,dict):deg=0
            return self.node_ins+(deg*self.edge_ins)+edge_weight*self.edge_weight_coef
        else:
            return sys.maxsize