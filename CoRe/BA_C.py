# Run gel modeling code
import sys
from optparse import OptionParser
import numpy as np
import math

class BA(object):
    def __init__(self):
        pass

    def set_response(self,_response):
        self.response_a = _response

    def get_CC(self):
        response = self.response_a

        rows, cols = response.shape[0], response.shape[1]

        # Assume uniform distribution as initial condition
        probs = np.ones(shape=(1,rows))/float(rows)

        # Tolerance for convergence
        errtol = 1e-4

        # Initial value for error to start the iterativr loop
        err = 1

        while err>errtol:
            c = np.zeros(shape=(rows,))

            v1 = np.matmul(probs,response)[0,:]

            for j in range(0,rows):
                for k in range(0,cols):
                    if response[j,k]>0.0:
                        c[j] += response[j,k]*math.log(response[j,k]/v1[k])

                c[j] = math.exp(c[j])

            mean_c = np.dot(probs,c)[0]

            I_L = math.log(mean_c)
            I_U = math.log(np.max(c))

            err = abs(I_U - I_L)

            if err>errtol:
                for j in range(0,rows):
                    probs[0,j] *= c[j]/mean_c
            else:
                C = I_L

        C *= 1.0/math.log(2)

        return C, err, probs
