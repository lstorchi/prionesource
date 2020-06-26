from numpy import diag, dot, float64, outer, sign, transpose, zeros
from numpy.linalg import det, svd, inv

import numpy

##############################################################################

def return_rotation_matrix (coord_from=None, coord_to=None):
  # Initialise the covariance matrix A.
  A = zeros((3, 3), float64)

  centroid_from = coord_from.sum(0) / coord_from.shape[0]
  centroid_to = coord_to.sum(0) / coord_to.shape[0]

  # Loop over the atoms.
  for i in range(coord_from.shape[0]):
    # The positions shifted to the origin.
    orig_from = coord_from[i] - centroid_from
    orig_to = coord_to[i] - centroid_to
                                
    # The outer product.
    A += outer(orig_from, orig_to)
                                          
  U, S, V = svd(A)
                                             
  # The handedness of the covariance matrix.
  d = sign(det(A))
  D = diag([1, 1, d])
                                                        
  R = dot(transpose(V), dot(D, transpose(U)))

  print("Verify rotation matrix ")
  print("det: (should be +/-1)" , det(R))
  print("transpose matrix should be equal to inverse matrix") 
  print(transpose(R))
  print(" ")
  print(inv(R))
            
  return R, centroid_from, centroid_to

##############################################################################

def rmsd(V, W):
  D = len(V[0])
  N = len(V)
  rmsd = 0.0
  for v, w in zip(V, W):
    rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
   
  return numpy.sqrt(rmsd/N)

##############################################################################

def centroid(X):
  C = sum(X)/len(X)

  return C
