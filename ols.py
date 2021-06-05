import numpy as np
import scipy.stats as st

"""
Adapted from https://github.com/sbird/spb_common/blob/master/leastsq.py

Module for computing the 5 least squares regression methods detailed in
        Linear regression in astronomy.
        Isobe, Takashi; Feigelson, Eric D.; Akritas, Michael G.; Babu, Gutti Jogesh
        Astrophysical Journal, Part 1 (ISSN 0004-637X), vol. 364, Nov. 20, 1990, p. 104-113
        http://adsabs.harvard.edu/abs/1990ApJ...364..104I
These methods are appropriate when the intrinsic scatter in the data is much larger
than the error on each data point.
"""

def leastsq(X,y, method=5):
    """
       Compute the least squares fit to y = beta x + alpha,
       using one of the 5 methods outlined in
       http://adsabs.harvard.edu/abs/1990ApJ...364..104I
       Method 1 minimises distance from Y given X (ie, the standard least squares fit)
       Method 2 minimises distance from X given Y
       Method 3 (recommended) is the OLS bisector, which gives a line bisecting the above two.
       Method 4 (Orthogonal regression) minimises perpendicular distance from the line to points
       Method 5 is the geometric mean of the slopes from methods 1 and 2.
       Method 6 is the Theil-Sen estimator: the median of the pairwise slopes.
       (See Akritas 95,  http://www.tandfonline.com/doi/abs/10.1080/01621459.1995.10476499)
       Returns:
              (alpha, beta, bvar), the intercept slope and variance of the slope
    """
    x = X[:,0]
    n = x.shape[0]
    k = 2

    #Define some sums
    xbar = np.mean(x) # Eq. 2a
    ybar = np.mean(y) # Eq. 2b
    xdif = x-xbar # Eq. 5
    ydif = y-ybar # Eq. 6
    sxx = np.sum(xdif**2) # Eq. 3a
    syy = np.sum(ydif**2) # Eq. 3b
    sxy = np.sum(ydif*xdif) # Eq. 4

    #Check for zeros
    if sxx == 0 or syy == 0 or sxy == 0:
        raise ValueError("Least Squares ill-defined")
    if method > 6 or method < 1:
        raise ValueError("Method not recognised")

    #These formulas are taken from Table 1 of Isobe et al, page 3
    #Minimise distance from Y given X
    beta1 = sxy/sxx # Table 1
    #Variance of b1
    bvar1 = np.sum(xdif**2*(ydif-beta1*xdif)**2)/sxx**2 # Table 1
    #Minimise distance from X given Y
    beta2 = syy/sxy # Table 1
    #Variance of b2
    bvar2 = np.sum(ydif**2*(ydif-beta2*xdif)**2)/sxy**2 # Table 1
    #Covariance of b1 and b2
    covb12 = np.sum(xdif*ydif*(ydif-beta2*xdif)*(ydif-beta1*xdif))/(beta1*sxx**2) # Note Table 1


    if method == 1:
        beta = beta1
        bvar = np.sqrt(bvar1)
        avar = np.sqrt(np.sum((ydif-beta*xdif-n*xbar*xdif*(ydif-beta*xdif)/sxx)**2)/n**2) # Eq. 9,11
    if method == 2:
        beta = beta2
        bvar = np.sqrt(bvar2)
        avar = np.sqrt(np.sum((ydif-beta*xdif-n*xbar*ydif*(ydif-beta*xdif)/sxy)**2)/n**2) # Eq. 9,15
    if method == 3:
        #OLS bisector: line that bisects the above two.
        beta1p1 = 1+beta1**2
        beta2p1 = 1+beta2**2
        beta = (beta1*beta2 - 1 + np.sqrt(beta1p1*beta2p1))/(beta1+beta2)
        #Variance
        prefac = beta**2 / ( (beta1 + beta2)**2 * beta1p1 * beta2p1)
        var = beta2p1**2 * bvar1 + 2 * beta1p1 * beta2p1 * covb12 + beta1p1**2 * bvar2
        bvar = np.sqrt(prefac*var)

        gamma1 = beta/((beta1+beta2)*np.sqrt(beta1p1*beta2p1)) # Eq. 20
        gamma13 = gamma1*beta2p1 # Eq. 12
        gamma23 = gamma1*beta1p1 # Eq. 17
        avar = np.sqrt(np.sum((ydif-beta*xdif-n*xbar*(gamma13*xdif*(ydif-beta1*xdif)/sxx + \
                                                      gamma23*ydif*(ydif-beta2*xdif)/sxy))**2)/n**2) # Eq. 9,12,17,20
    """
    TO DO: Implement avar for the methods 4 and 5.
    Check reference for the method 6.
    """
    if method == 4:
        #Orthogonal: minimise perpendicular distance from line to points
        beta = 0.5*((beta2-1./beta1)+np.sign(sxy)*np.sqrt(4+(beta2-1./beta1)**2))
        prefac = beta**2 / (4*beta1**2 + (beta1*beta2 - 1)**2)
        bvar = prefac * ( bvar1/beta1**2 + 2*covb12 + beta1**2*bvar2 )
        avar = 0

    if method == 5:
        #Reduced major axis:
        beta = np.sign(sxy)*np.sqrt(beta1*beta2)
        bvar = 0.25 * (beta2/beta1 * bvar1 + 2*covb12 + beta1/beta2 * bvar2)
        avar = 0

    if method == 6:
        #Theil-Sen estimator for uncensored data: the median of the slopes.
        yy = np.subtract.outer(y,y)
        xx = np.subtract.outer(x,x)
        ind = np.where(xx != 0)
        beta = np.median(yy[ind]/xx[ind])
        #Can't find a formula for the variance
        bvar = 0
        avar = 0

    #The intercept
    alpha = ybar - beta*xbar
    yhat = alpha + beta*x
    rms = np.sqrt(sum((y-yhat)**2)/(n-k)) # in units of Y

    return (alpha, avar, beta, bvar, rms)

def pearson(x,y,alpha, beta, method=3):
    """Find the Pearson correlation coefficient between the fit and the data.
        if method == 1 return the Pearson r of y and the fit to y
        if method == 2, the same but with x and y reverse
        if method == 3 the geometric mean of the above
    """
    #Vector of expected y from fit
    fity = beta*x + alpha
    #Vector of expected x from fit
    fitx = (y - alpha) / beta
    #Scatter from y axis: method 1 minimises this.
    (pry,_) = st.pearsonr(y,fity)
    if method == 1:
        return pry
    (prx,_) = st.pearsonr(x,fitx)
    if method == 2:
        return prx
    return np.sqrt(pry*prx)

def kstest(x,y,alpha, beta):
    """Find the K-S test probability that the fit and the data were from the same distribution"""
    #Vector of expected y from fit
    fity = beta*x + alpha
    #Vector of expected x from fit
    fitx = (y - alpha) / beta
    (D1, p1) = st.ks_2samp(y,fity)
    (D2, p2) = st.ks_2samp(x,fitx)
    return (np.sqrt(D1*D2),np.sqrt(p1*p2))
