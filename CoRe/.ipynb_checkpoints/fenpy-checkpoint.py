#!/usr/bin/env python
# coding: utf-8
# %%
"""Contains functions for calculating the propagation of entropy and state relative entropies for graphs
"""
import numpy as np
import math

def state_entropy(state,sources):
    H = 0.0

    for i in range(0,state.shape[0]):
        if int(i/2) not in sources:
            if state[i]>0.0:
                H += state[i]*math.log(state[i]/0.5)

    H *= 1.0/math.log(2)

    return H

def get_final_state(Tmat,X,sources):
    tol = 1.0

    while tol>1e-6:
        X_old = np.array(X, copy=True)
        X = np.matmul(Tmat,X)

        for s in sources:
            X[2*s:2*(s+1),:] = X_old[2*s:2*(s+1),:]

        tol = np.linalg.norm(X-X_old,2)

        #print(np.min(X),np.max(X))

    return X
