import numpy as np

def lhs_ellipse(N_prop, axes):
    """
    This function samples approximatly N_prop data points in a LHS method
    from an ellipse with the ellipsoid axes 'axes'.
    This is done computation and memory efficient.
    """
    # dimensionality of the ellipsoid
    dim = len(axes)

    # expected volume of a n-sphere
    V_exp = np.pi**(dim/2.0) / np.math.gamma(dim/2 + 1)*axes.prod()

    # the expected number of samples within a cube
    N_req = int(N_prop *  axes.prod() * 2**dim / V_exp)

    # initialize the loop quantities
    mask = np.ones(N_req)
    result = np.ones((N_req,1))

    # loop over each dimension and remove for each new dimension excluded points.
    # think of a 2 dim projection of a 3 dim sphere. On the edges of the square points can be removed.
    # this is done here iterativly from first dimension up to dimension 'dim'
    for i in range(dim):
        print(i)
        distance = np.zeros(N_req)

        if i!=0:
            result = np.hstack([result, np.empty([N_req,1])])

        # propose values of new dimension
        proposal = np.linspace(-axes[i],axes[i],N_req,endpoint=False) + np.random.rand(N_req)*2*axes[i] / N_req
        np.random.shuffle(proposal)
        result[:,i] = proposal

        # check whether some points lay outside the ellipsoid
        for j in range(i+1):
            distance+=result[:,j]**2/axes[j]**2
        mask = distance < 1
        result =result[mask]

        # update the size of points which are still in the game. It will comverge to 'N_prop' after all dimensions.
        N_req=len(result)

    return(result)