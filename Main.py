#!use/bin/env python
# -*- coding:utf-8 -*-

'''
Created on 20220516

@author: think
'''
from inverse.ARAP import ARAP
from inverse.ARAP import cellNewSolution,rotNewSolution
import numpy as np
from math import pi
from copy import deepcopy

if __name__ == '__main__':
    a = ARAP()
    p1 = [1957,133,0]
    p3 = [120,1950,0]
    p4 = [1750,1750,0]
    pm = [1023,1023,0]
    p1N = [1956.02175,133.02264,0]
    p3N = [118.977535,1949.97811,0]
    p4N = [1748.98238,1750.01762,0]
    pmN = deepcopy(pm)
    cell = a.cellCreate(p1,p3,p4,pm)
    premove = a.cellCreate(p1N,p3N,p4N)
    #print(a.calRot(cell,cellNew))
    
    '''
    pmove = a.setPMove(cell,remove=[1])
    omega = a.setCellOmega(cell.shape[0])
    edge = a.setEdge(cell)
    r = np.tile(np.eye(3),(len(cell),1,1))
    #L = a.setPLMatrix(omega)
    removeIndex = [3]
    L,b = a.setSolutionMatrix(omega, r, edge, removeindex=removeIndex, premove=premove)
    #L,b = a.setSolutionMatrix(omega, r, edge)
    print('the rank of L matrix is',np.linalg.matrix_rank(L))
    cellnew = np.insert(cellNewSolution(L, b), removeIndex, premove, axis=0)
    print(cellnew)
    rotnew = rotNewSolution(cell, cellnew)
    print(rotnew)
    '''
    removeIndex = [0,1,2]
    r = np.tile(np.eye(3),(len(cell),1,1))
    cellguess = a.cellCreate(p1N,p3N,p4N,pmN)
    rotnew, cellnewinternew, count = a.calInteration(cell, cellguess, removeindex = removeIndex, premove=premove)
    disp = cellnewinternew - cell
    print('interation num is:', count)
    print('rot in points is:')
    for ele in rotnew:
        print('{:.3f}\n'.format(ele[1,0]*180*3600/pi))
    print('disp in points is:')
    for ele in disp:
        print('\n{:.5f}  {:.5f}  {:.5f}'.format(*ele))
    
    
    
    