#!usr/bin/env python
# -*- coding:utf-8 -*-

'''
Created on 20220514

@author: think
'''

import numpy as np
from _testbuffer import ndarray

class ARAP:
    def __init__(self):
        self.R = 0
    def calInteration(self, cellinit, cellguess, **kwargs):
        cellfore = cellinit#p0
        cellnewinternew = cellguess#p0'guessed
        count = 0
        countThre = 1e2
        rotthreshold = 1e-6
        cellthreshold = 1e-6
        if kwargs:
            removeIndex = kwargs['removeindex']
            premove = kwargs['premove']
        omega = self.setCellOmega(cellfore.shape[0])
        edge = self.setEdge(cellfore)
        
        ##############################
        #calculate rot new  r0
        ##############################
        rotnew = rotNewSolution(cellfore, cellnewinternew)#r0(p0,p0')
        
        
        ##############################
        #calculate cell new  p1'
        ##############################
        cellnewinterfore = cellnewinternew#p0'
        L, b = self.setSolutionMatrix(omega, rotnew, edge, removeindex=removeIndex, premove=premove)
        cellnewinternew = cellNewSolution(L, b)
        cellnewinternew = self.cellInsert(cellnewinternew, removeIndex, premove)#p1'(p0,r0)
        
        
        ##############################
        #calculate rot new  r1
        ##############################
        rotfore = rotnew#r0
        rotnew = rotNewSolution(cellfore, cellnewinternew)#r1(p0,p1')
        
        
        
        ##############################
        #interation
        ##############################
        rotdelte = self.calRotDelte(rotfore, rotnew)#r1-r0
        celldelte = self.calCellDelte(cellnewinterfore, cellnewinternew)#p1'-p0'
        
        
        while np.amax(rotdelte)>rotthreshold or np.amax(celldelte)>cellthreshold:
            ##############################
            #calculate cell new
            ##############################
            cellnewinterfore = cellnewinternew
            L, b = self.setSolutionMatrix(omega, rotnew, edge, removeindex=removeIndex, premove=premove)
            cellnewinternew = cellNewSolution(L, b)
            cellnewinternew = self.cellInsert(cellnewinternew, removeIndex, premove)
            
            
        
            ##############################
            #calculate rot new
            ##############################
            rotfore = rotnew
            rotnew = rotNewSolution(cellfore, cellnewinternew)
        
        
        
            ##############################
            #interation
            ##############################
            rotdelte = self.calRotDelte(rotfore, rotnew)
            celldelte = self.calCellDelte(cellnewinterfore, cellnewinternew)
            
            count += 1
            if count > countThre:
                break
        
        return rotnew, cellnewinternew, count
    
    def calRot(self, omegaEdge, edge, edgeNew):#calculate rot of one certain vertex==
        '''
        omegaEdge = [omega1,omega2,omega3] = self.setOmega()
        edge = [edge1,edge2,edge3] = self.setEdge(cell)
        edgeNew = [edgeNew1,edgeNew2,edgeNew3] = self.setEdge(cellNew)
        '''
        
        '''
        ###################
        omegaEdge:1*n
        edge:1*n
        ###################
        '''
        SMatrix = np.dot(omegaEdge*edge.T, edgeNew)
        [u,s,v] = np.linalg.svd(SMatrix)
        rot = np.dot(v,u.T)
        return rot
    def setEdge(self,cell):
        '''
        ################
        edgeMatrix n*n
        #################
        '''
        edgeMatrix = []
        for i in range(len(cell)):
            edgeMatrix.append([])
            for j in range(len(cell)):
                edgeMatrix[i].append(cell[i]-cell[j])
        return np.array(edgeMatrix)
        '''
        if cell.shape == (4,3):
            return np.array([cell[0][:]-cell[1][:] , cell[0][:]-cell[2][:] , cell[0][:]-cell[3][:]]).T
        else:
            print('error p shape')
        '''
    def setOmega(self):
        w1 = np.array([1,1,1])
        w2 = np.array([1,1,1])
        w3 = np.array([1,1,1])
        return np.array([w1,w2,w3]).T
        #设置移动点的系数
    
    def setCellOmega(self, cellnum:int)->ndarray:
        '''
        ###################
        [[0,1,1,1],
         [1,0,1,1],
         [1,1,0,1],
         [1,1,1,0]] n*n
        ##################
        '''
        cellOmega = np.ones([cellnum,cellnum])
        for i in range(cellnum):
            cellOmega[i,i] = 0
        return np.array(cellOmega)
    def cellCreate(self,*p):
        cell = np.array(list(p))
        return cell
    def setPMove(self, cell, **premove)->ndarray:
        pmove = np.array(list(cell))
        premove = premove['remove']
        pmove = np.delete(pmove,premove,axis=0)
        return pmove
    def setPLMatrix(self,cellomega,**kwargs):
        L = -cellomega
        for i in range(len(cellomega)):
            L[i,i] = np.sum(cellomega[i,:]) - cellomega[i,i]
        if kwargs:
            removeIndex = kwargs['removeindex']
            Lremove = np.delete(L,removeIndex,axis=0)
            Lremove = np.delete(Lremove,removeIndex,axis=1)
            return L,Lremove
        return L
    def setSolutionMatrix(self, cellomega, r, edge, **kwargs):
        rMatrix = []
        for i in range(len(r)):
            rMatrix.append([])
            for j in range(len(r)):
                rMatrix[i].append(r[i]+r[j])
        rMatrix = np.array(rMatrix)
        bMatrix = []
        for i in range(len(r)):
            bMatrix.append([])
            for j in range(len(r)):
                bMatrix[i].append(np.dot((cellomega[i,j]*rMatrix[i,j]),edge[i,j]).tolist())
        pb = 0.5*np.sum(np.array(bMatrix),axis = 1)
        if kwargs:
            L,Lremove = self.setPLMatrix(cellomega, removeindex=kwargs['removeindex'])
            removeIndex = kwargs['removeindex']
            premove = kwargs['premove']
            boolIndex = np.zeros(len(r), dtype = bool)
            boolIndex[removeIndex] = True
            pb = np.delete(pb,removeIndex,axis=0)
            for i in range(len(removeIndex)):
                temp = L[:,removeIndex[i]]
                pb -= temp[~boolIndex].reshape(len(Lremove),)*premove[i]
        else:
            print('no control point! lack rank matrix')
            Lremove = self.setPLMatrix(cellomega)
        return Lremove,pb
    '''
    def rotNewSolution(self, cell,cellnew):
        omega = self.setCellOmega(self, cell)
        edge = self.setEdge(self, cell)
        edgeNew = self.setEdge(self, cellnew)
        rmatrix = []
        for i in range(len(cell)):
            rmatrix.append(self.calRot(self, omega, edge, edgeNew))
    '''
    def calRotDelte(self, rotfore, rotnew):
        return self.calDelta(rotfore, rotnew)
    def calCellDelte(self, cellfore, cellnew):
        return self.calDelta(cellfore, cellnew)
    def calDelta(self, before, new):
        delta = []
        for i in range(len(before)):  
            delta.append(np.sum(np.power((before[i]-new[i]),2).tolist()))
        return np.array(delta)
    def cellInsert(self, cell, removeindex, premove):
        for i in range(len(removeindex)):
            cell = np.insert(cell, removeindex[i], premove[i], axis=0)
        return cell
def cellNewSolution(pl,pb):
    cellnew = np.linalg.solve(pl,pb)
    return cellnew
def rotNewSolution(cell,cellnew):
    omega = ARAP.setCellOmega(ARAP,cell.shape[0])
    edge = ARAP.setEdge(ARAP,cell)
    edgeNew = ARAP.setEdge(ARAP,cellnew)
    rmatrix = []
    for i in range(len(cell)):
        rmatrix.append(ARAP.calRot(ARAP,omega[i], edge[i], edgeNew[i]))
    return np.array(rmatrix)
        