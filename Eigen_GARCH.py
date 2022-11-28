import numpy as np
import pandas as pd

def eigenLikelihood(parameters, x, lambdas, out=None):
    """ Function for the Eigen GARCH introduced in Hetland et al(2019), with three return series:
    Parameters
    ----------
    paramters: list : the 24 parameters to be estimated in the model
    x: numpy array : a 3xN array of the three return series
    lambdas: numpy array: a 3xN array of the eigenvalues
            
    Returns
    =======
    llh: the loglikelihood of the estimation
    """
    # Dimensions of the return series
    p,T=x.shape 
    
    # Unpack initial values for the parameters in the model
    phi1, phi2, phi3, w1, w2, w3, a11, a12, a13, a21, a22, a23, a31, a32, a33, b11, b12, b13, b21, b22,b23, b31, b32,b33 = parameters
    
    # Use parameters to define the variables in the model
    W = np.array(([w1], [w2], [w3]))
    A = np.array(([a11, a12,a13],
                  [a21, a22,a23],
                  [a31, a32,a33]))
    B = np.array(([b11, b12,b13],
                  [b21, b22,b23],
                  [b31, b32,b33]))
    V = np.array(([np.cos(phi1), np.sin(phi1),0],[-np.sin(phi1), np.cos(phi1),0 ], [0, 0,1 ])
                )@ np.array(([np.cos(phi2), 0, np.sin(phi2)], [0, 1,0 ], [-np.sin(phi2), 0, np.cos(phi2)])
                           )@ np.array(( [1, 0,0 ], [0,np.cos(phi3), np.sin(phi3)], [0,-np.sin(phi3), np.cos(phi3)]))
    Xtilde = V.transpose() @ x

    # Conditional eigenvalue dynamics    
    lambdas = (np.ones((T,p))).transpose()

    for t in range(1,T):
        lambdas[:,t:t+1] = W + A @ np.multiply(Xtilde[:,t-1:t],Xtilde[:,t-1:t]) + B @ lambdas[:,t-1:t]

    # Loglikelihoods contributions for each return series
    lls1 = np.log(lambdas[0,:]) + np.multiply(Xtilde[0,:],Xtilde[0,:]) / lambdas[0,:] 
    lls2 = np.log(lambdas[1,:]) + np.multiply(Xtilde[1,:],Xtilde[1,:]) / lambdas[1,:] 
    lls3 = np.log(lambdas[2,:]) + np.multiply(Xtilde[2,:],Xtilde[2,:]) / lambdas[2,:] 
    lls = (lls1 + lls2 + lls3)*p*np.log(2*np.pi)
    
    # Loglikelihood
    ll = 0.5*(np.sum(lls1 + lls2+lls3) + T*p*np.log(2*np.pi))
    # Return
    return ll

def eigenresiduals(parameters, x, lambdas, out=None):
    """ Function for the Eigen GARCH introduced in Hetland et al(2019), with three return series:
    Parameters
    ----------
    paramters: list : the 24 parameters to be estimated in the model
    x: numpy array : a 3xN array of the three return series
    lambdas: numpy array: a 3xN array of the eigenvalues
            
    Returns
    =======
    llh: the loglikelihood of the estimation
    """
    # Dimensions of the return series
    p,T=x.shape 
    
    # Unpack initial values for the parameters in the model
    phi1, phi2, phi3, w1, w2, w3, a11, a12, a13, a21, a22, a23, a31, a32, a33, b11, b12, b13, b21, b22,b23, b31, b32,b33 = parameters
    
    # Use parameters to define the variables in the model
    W = np.array(([w1], [w2], [w3]))
    A = np.array(([a11, a12,a13],
                  [a21, a22,a23],
                  [a31, a32,a33]))
    B = np.array(([b11, b12,b13],
                  [b21, b22,b23],
                  [b31, b32,b33]))
    V = np.array(([np.cos(phi1), np.sin(phi1),0],[-np.sin(phi1), np.cos(phi1),0 ], [0, 0,1 ])
                )@ np.array(([np.cos(phi2), 0, np.sin(phi2)], [0, 1,0 ], [-np.sin(phi2), 0, np.cos(phi2)])
                           )@ np.array(( [1, 0,0 ], [0,np.cos(phi3), np.sin(phi3)], [0,-np.sin(phi3), np.cos(phi3)]))
    Xtilde = V.transpose() @ x

    # Conditional eigenvalue dynamics    
    lambdas = (np.ones((T,p))).transpose()
    residuals = (np.ones((T,p))).transpose()
    
    for t in range(1,T):
        lambdas[:,t:t+1] = W + A @ np.multiply(Xtilde[:,t-1:t],Xtilde[:,t-1:t]) + B @ lambdas[:,t-1:t]
        residuals[:,t:t+1] = np.multiply(np.sqrt(np.diag(lambdas[:,t:t+1])),Xtilde[:,t-1:t])
        
    # Residuals
    
    # Return
    return residuals