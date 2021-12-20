import numpy as np
import pandas as pd
from scipy.special import gammaln, psi, factorial
from scipy.optimize import fmin_l_bfgs_b as optim
from scipy.stats import nbinom
from sklearn.mixture import GaussianMixture as GMM
import time
import multiprocessing as mp
import sys

# X is a numpy array representing the data
# initial params is a numpy array representing the initial values of size and prob parameters
def nbfit(X):
    infinitesimal = np.finfo(np.float).eps
    def log_likelihood(params, *args):
        r, p = params
        if r == 0:
            r = infinitesimal
        X = args[0]
        N = X.size
        #MLE estimate based on the formula on Wikipedia:
        # http://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation (correct formula for Wikipedia terms for f(k;r,p) r-number of failures, k-number of successes).
        # The definition of the negative binomial distribution can be extended to the case where the parameter r can take on a positive real value.
        # The maximum likelihood estimator only exists for samples for which the sample variance is larger than the sample mean.
        result = np.sum(gammaln(X + r)) \
            - np.sum(np.log(factorial(X))) \
            - N*(gammaln(r)) \
            + N * r * np.log((1-p if p < 1 else infinitesimal)) \
            + np.sum(X * np.log(p))
            # Both Wikipedia and np.negative_binomial have k and r swapped.
            # + N*r*np.log(p) \
            # + np.sum(X*np.log(1-(p if p < 1 else 1-infinitesimal)))
        return -result

    def log_likelihood_deriv(params, *args):
        r, p = params
        if r == 0:
            r = infinitesimal
        X = args[0]
        N = X.size
        #pderiv = (N*r)/p - np.sum(X)/(1-(p if p < 1 else 1-infinitesimal))
        #rderiv = np.sum(psi(X + r)) - N*psi(r) + N*np.log(p)
        pderiv = (N*r)/(p if p < 1 else 1-infinitesimal) - np.sum(X)/(1-p if p < 1 else infinitesimal)
        rderiv = np.sum(psi(X + r)) - N*psi(r) + N*np.log(p if p < 1 else 1-infinitesimal)
        return np.array([-rderiv, -pderiv])

    #reasonable initial values (from fitdistr function in R)
    m = np.mean(X) #Estimation of number of failures.
    v = np.var(X)
    size = (m**2)/(v-m) if v > m else 10
    #convert size,mu to size,probability of success.
    p0 = size / ((size+m) if size+m != 0 else 1)
    r0 = size
    initial_params = np.array([r0, p0])
    #Minimize a function func using the L-BFGS-B algorithm.
    # x0: Initial guess.
    # args: sequence (optional) Arguments to pass to func and fprim
    # approx_gradbool: Whether to approximate the gradient numerically (in which case func returns only the function value).
    # bounds: (min, max) pairs for each element in x, defining the bounds on that parameter. Use None or +-inf for one of min or max when there is no bound in that direction.
    bounds = [(infinitesimal, None), (infinitesimal, 1)] #Bounds for n and p.
    optimres = optim(log_likelihood, x0=initial_params, args=(X + infinitesimal,), fprime=log_likelihood_deriv, bounds=bounds) #, approx_grad=1
    params = optimres[0]
    return {'size': params[0], 'prob': params[1]}


#Calculate UDP for each pair of probe,sample. UDP is the probability of the higher Normal distribution.
def calc_udp_nbm(data):
    data = data.apply(lambda row: row.fillna(row.mean()), axis=1)
    my_udp = np.empty((0, len(data.columns)))
    for index, row in data.iterrows():
        row = row.astype(int).values
        #Remove outliers
        if(max(row)<5 or np.count_nonzero(row)<10): #(row == 0).all() or np.sum(row)<200
            my_udp = np.append(my_udp, [np.zeros(len(data.columns))], axis=0)
            continue
        row[np.argsort(row)[-3]] = row[np.argsort(row)[-2]] = row[np.argsort(row)[-1]] = row[np.argsort(row)[-4]]
        res = nbfit(row)
        size = res['size']
        prob = res['prob']
        row_probs = nbinom.pmf(row, size, prob)
        my_udp = np.append(my_udp, [row_probs], axis=0)
    return(pd.DataFrame(data=my_udp, index=data.index, columns=data.columns))


def calc_udp_gmm(data):
    data = data.apply(lambda row: row.fillna(row.mean()), axis=1)
    my_udp = np.empty((0, len(data.columns)))
    for index, row in data.iterrows():
        row = row.values.reshape(-1,1)
        gmm = GMM(n_components=2, covariance_type='full', random_state=0).fit(row)
        pred = gmm.predict_proba(row)[:,0]
        my_udp = np.append(my_udp, [pred], axis=0)
    return(pd.DataFrame(data=my_udp, index=data.index, columns=data.columns))


#Run calc_udp on parallel.
def calc_udp_multi_process(rma_filename, file_type = 0):
    print(time.ctime(), 'Calculate UDP...')
    if (file_type == 0):
        udp_func = calc_udp_gmm
    else:
        udp_func = calc_udp_nbm
    udp = pd.DataFrame()
    rma_file = pd.read_csv(rma_filename, index_col=0, chunksize=700)
    pool = mp.Pool(6)  # use 6 processes
    results = [pool.apply_async(udp_func, args=(x,)) for x in rma_file]
    for p in results:
        udp = udp.append(p.get())  # f.get(timeout=100)
        print('.', end="")
        sys.stdout.flush()
    udp.to_csv('data/output_udp.csv')


#Testing of the negative binom fitting function.
if __name__ == '__main__':
    """
    Samples are drawn from a negative binomial distribution with specified parameters,
    n successes and p probability of success where n is an integer > 0 and p is in the
    interval [0, 1]. Drawn samples where each sample is equal to N, the number of
    failures that occurred before a total of n successes was reached. N+n is the number
    of trials. The distribution P(N) gives the probability of N failures given n
    successes, with a success on the last trial.
    #This is not exactly the same as the Wikipedia definition.
    # Series of N1,N2,  (failures until n successes appear).
    # P(X=N) = (N+n-1)over(N) * p^n * (1-p)^N       (different from the Wikipedia formula by the base of the binomial).
    #testset = np.random.negative_binomial(float(size), float(prob), 1000)
    https://en.wikipedia.org/wiki/Negative_binomial_distribution#Probability_mass_function
    """


    size = 10 #r-successes until N failures happen.
    prob = 0.3 #probability of success.
    testset = nbinom.rvs(size, prob, size=1000)
    res = nbfit(testset)
    _size = res['size']
    _prob = res['prob']
    row_probs = nbinom.pmf(testset, size, prob)
    print(f'Size: {_size}, prob: {_prob}')
    #print(f'Probabilities: {row_prob}')
    #print(f'Error size: {error_size}, error prob: {error_prob}')


    """
    use_udp = input('Use pre-calculated UDP file? (y/n):')
    if(use_udp == 'y'):
        udp_filename = input('Enter UDP filename:')
        udp = pd.read_csv(udp_filename, index_col=0)
    else:
        rma_filename = input('Enter RMA filename:') #'src/Pathologist/Input/Lung/GSE29013_RMA.txt'
        udp = calc_udp_multi_process(rma_filename)
     """

    """
    udp = calc_udp_multi_process('ACC.rnaseq.txt')
    results, _ = calc_activity_consistency(udp)
    results.to_csv('data/output_activity.csv', index=False)
    print(time.ctime(), 'Done.')
    """