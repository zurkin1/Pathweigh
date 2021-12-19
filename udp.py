# Fitting various statistical distributions to either microarray or RNAseq data.
# In general continous distributions have a fit function, and so we use it to fit the data.
# From documentation, it is 'maximum likelihood estimation of distribution parameters, including location
# and scale'.
# Discrete distributions don't have a fit function hence we use either a closed form formula
# For maximum likelihood or an iterative optimization process.
import numpy as np
import pandas as pd
from scipy.special import gammaln, psi, factorial
from scipy.optimize import fmin_l_bfgs_b as optim
from scipy.stats import nbinom, norm, poisson, gennorm
from sklearn.mixture import GaussianMixture as GMM
from numpy import inf
import gc
import time
import multiprocessing as mp
import sys


infinitesimal = np.finfo(np.float).eps
relative_path = './Pathweigh/'


def log_likelihood_nb(params, *args):
    r, p = params
    if r == 0:
        r = infinitesimal
    X = args[0]
    N = X.size
    # MLE estimate based on the formula on Wikipedia:
    # http://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation (correct formula for Wikipedia terms for f(k;r,p) r-number of failures, k-number of successes).
    # The definition of the negative binomial distribution can be extended to the case where the parameter r can take on a positive real value.
    # The maximum likelihood estimator only exists for samples for which the sample variance is larger than the sample mean.
    result = np.sum(gammaln(X + r)) \
        - np.sum(np.log(factorial(X))) \
        - N * (gammaln(r)) \
        + N * r * np.log((1 - p if p < 1 else infinitesimal)) \
        + np.sum(X * np.log(p))
    # Both Wikipedia and np.negative_binomial have k and r swapped.
    # + N*r*np.log(p) \
    # + np.sum(X*np.log(1-(p if p < 1 else 1-infinitesimal)))
    return -result


def log_likelihood_deriv_nb(params, *args):
    r, p = params
    if r == 0:
        r = infinitesimal
    X = args[0]
    N = X.size
    # pderiv = (N*r)/p - np.sum(X)/(1-(p if p < 1 else 1-infinitesimal))
    # rderiv = np.sum(psi(X + r)) - N*psi(r) + N*np.log(p)
    pderiv = (N * r) / (p if p < 1 else 1 - infinitesimal) - np.sum(X) / (1 - p if p < 1 else infinitesimal)
    rderiv = np.sum(psi(X + r)) - N * psi(r) + N * np.log(p if p < 1 else 1 - infinitesimal)
    return np.array([-rderiv, -pderiv])

# X is a numpy array representing the data
# initial params is a numpy array representing the initial values of size and prob parameters


def nbfit(X):
    # reasonable initial values (from fitdistr function in R)
    m = np.mean(X)  # Estimation of number of failures.
    v = np.var(X)
    size = (m ** 2) / (v - m) if v > m else 10
    # convert size,mu to size,probability of success.
    p0 = size / ((size + m) if size + m != 0 else 1)
    r0 = size
    initial_params = np.array([r0, p0])
    # Minimize a function func using the L-BFGS-B algorithm.
    # x0: Initial guess.
    # args: sequence (optional) Arguments to pass to func and fprim
    # approx_gradbool: Whether to approximate the gradient numerically (in which case func returns only the function value).
    # bounds: (min, max) pairs for each element in x, defining the bounds on that parameter. Use None or +-inf for one of min or max when there is no bound in that direction.
    bounds = [(infinitesimal, None), (infinitesimal, 1)]  # Bounds for n and p.
    optimres = optim(log_likelihood_nb, x0=initial_params, args=(X + infinitesimal,), fprime=log_likelihood_deriv_nb,
                     bounds=bounds)  # , approx_grad=1
    params = optimres[0]
    # f_value = optimres[1] #Maximm likelihood reached.
    # d['warnflag'] is 0 if converged, 1 if too many function evaluations or too many iterations,
    # 2 if stopped for another reason, given in d['task'].
    # warn = optimres[2] #Dictinary with convergance resutls.
    return {'size': params[0], 'prob': params[1]}  # , 'f_value': f_value, 'warn': warn, 'aic':aic}


# Calculate UDP for each pair of probe,sample. UDP is the probability of the higher Normal distribution.
def calc_udp_nbm(data, aic_test=False):
    my_udp = np.empty((0, len(data.columns)))
    aic = 0
    for index, row in data.iterrows():
        row = row.astype(int).values
        # Remove outliers
        if (max(row) < 5 or np.count_nonzero(row) < 10):  # (row == 0).all() or np.sum(row)<200
            my_udp = np.append(my_udp, [np.zeros(len(data.columns))], axis=0)
            continue
        row[np.argsort(row)[-3]] = row[np.argsort(row)[-2]] = row[np.argsort(row)[-1]] = row[np.argsort(row)[-4]]
        res = nbfit(row)
        size = res['size']
        prob = res['prob']
        row_probs = nbinom.pmf(row, size, prob)
        my_udp = np.append(my_udp, [row_probs], axis=0)
        if(aic_test):
            row[row == 0] = 1
            LogLik = (np.log(nbinom.pmf(row, size, prob) + infinitesimal)).sum()
            aic += (2 * 2 - 2 * LogLik) / 1000  # (-log_likelihood((size, prob), row))
    return (pd.DataFrame(data=my_udp, index=data.index, columns=data.columns)), aic


def calc_udp_gmm(data, aic_test=False):
    my_udp = np.empty((0, len(data.columns)))
    aic = 0
    for index, row in data.iterrows():
        row = row.values.reshape(-1, 1)
        gmm = GMM(n_components=2, covariance_type='full', random_state=0).fit(row)
        pred = gmm.predict_proba(row)[:, 0]
        my_udp = np.append(my_udp, [pred], axis=0)
        if(aic_test):
            aic += gmm.aic(row) / 1000
    return (pd.DataFrame(data=my_udp, index=data.index, columns=data.columns)), aic


# https://stackoverflow.com/questions/20011122/fitting-a-normal-distribution-to-1d-data
def calc_udp_norm(data, aic_test=False):
    my_udp = np.empty((0, len(data.columns)))
    aic = 0
    k = 2  # len(fitted_params)
    for index, row in data.iterrows():
        if (row == 0).all():
            my_udp = np.append(my_udp, [np.zeros(len(data.columns))], axis=0)
            continue
        row = row.values.reshape(-1, 1)
        mu, std = norm.fit(row)
        row_probs = norm.cdf(row, mu, std)[:, 0]
        my_udp = np.append(my_udp, [row_probs], axis=0)
        if(aic_test):
            logLik = np.sum(norm.logpdf(row, loc=mu, scale=std))
            aic += (2 * k - 2 * (logLik)) / 1000
    return (pd.DataFrame(data=my_udp, index=data.index, columns=data.columns)), aic


def calc_udp_poisson(data, aic_test=False):
    # define a likelihood function. https://www.statlect.com/fundamentals-of-statistics/Poisson-distribution-maximum-likelihood
    def loglikelihood_f(lmba, x, neg=1):
        # Using Stirling formula to avoid calculation of factorial.
        # logfactorial(n) = n*ln(n) - n
        n = x.size
        logfactorial = np.log(factorial(x))  # x*np.log(x) - x #
        logfactorial[logfactorial == -inf] = 0
        result =\
            - np.sum(logfactorial) \
            - n * lmba \
            + np.log(lmba) * np.sum(x)
        return neg * result

    my_udp = np.empty((0, len(data.columns)))
    aic = 0
    for index, row in data.iterrows():
        row = row.astype(int).values
        # Remove outliers
        if (max(row) < 5 or np.count_nonzero(row) < 10):  # (row == 0).all() or np.sum(row)<200
            my_udp = np.append(my_udp, [np.zeros(len(data.columns))], axis=0)
            continue
        row[np.argsort(row)[-3]] = row[np.argsort(row)[-2]] = row[np.argsort(row)[-1]] = row[np.argsort(row)[-4]]
        # res = poisson.fit(row) #No fit for discrete scipy.stats distribution.
        #result = []
        # for i in range(np.mean(row)-60, np.mean(row)+60):  # in fact (80, 120) should probably be enough
        #    _ = optim.fmin(likelihood_f, [i, 0.5, 0], args=(X, -1), full_output=True, disp=False)
        #    result.append((_[1], _[0]))
        #P2 = sorted(result, key=lambda x: x[0])[0][1]
        # Poissong maximum likelihood estimator has a closed form formula.
        lmba = np.mean(row)
        row_probs = poisson.pmf(row, lmba)
        my_udp = np.append(my_udp, [row_probs], axis=0)
        if(aic_test):
            row[row <= 0] = 1
            aic += (2 * 1 - 2 * loglikelihood_f(lmba, row)) / 1000
    return (pd.DataFrame(data=my_udp, index=data.index, columns=data.columns)), aic


def calc_udp_gennorm(data, aic_test=False):
    my_udp = np.empty((0, len(data.columns)))
    aic = 0
    k = 1  # len(fitted_params)
    for index, row in data.iterrows():
        if (max(row) < 5 or np.count_nonzero(row) < 10):
            my_udp = np.append(my_udp, [np.zeros(len(data.columns))], axis=0)
            continue
        row = row.values.reshape(-1, 1)
        beta, loc, scale = gennorm.fit(row)  # https://en.wikipedia.org/wiki/Generalized_normal_distribution. gennorm.numargs==1.
        row_probs = gennorm.pdf(row, beta, loc=loc, scale=scale)[:, 0]
        my_udp = np.append(my_udp, [row_probs], axis=0)
        if(aic_test):
            logLik = np.sum(gennorm.logpdf(row, beta, loc=loc, scale=scale))
            aic += (2 * k - 2 * (logLik)) / 1000
    return (pd.DataFrame(data=my_udp, index=data.index, columns=data.columns)), aic


# Read chunks of dataframe rows.
def chunker_rows(data, size):
    len_df = len(data)
    return [data.iloc[pos:min(pos + size, len_df)] for pos in range(0, len_df, size)]


# Run calc_udp on parallel.
def calc_udp_multi_process(mix, is_rnaseq):
    data = mix #pd.read_csv('input.csv', index_col=0)
    gc.collect()
    print(time.ctime(), f'Calculate UDP, is_rnaseq: {is_rnaseq}')
    data = data.apply(lambda row: row.fillna(row.mean()), axis=1)
    probes = pd.read_csv(relative_path + 'probelinksv2.txt', index_col=0)
    genes = pd.read_csv(relative_path + 'gene_list.csv')
    if (not is_rnaseq):
        data = data.loc[list(set(data.index) & set(probes.probe.values))]
        func = calc_udp_gmm
    else:
        genes['gene'] = genes.gene.map(str.lower)
        data = data.loc[list(set(data.index) & set(genes.gene.values))]
        func = calc_udp_norm #nbm
    df = pd.DataFrame()
    pool = mp.Pool()  # Use number of CPUs processes.
    results = [pool.apply_async(func, args=(x,)) for x in chunker_rows(data, 700)]
    for p in results:
        task_df, aic = p.get()
        df = df.append(task_df)  # f.get(timeout=100)
        print('.', end="")
        sys.stdout.flush()
    df.to_csv('output_udp.csv')
    print(time.ctime(), 'Done.')
    return df


if __name__ == '__main__':
    #data = pd.read_csv('./data/GSE29013_RMA.csv', index_col=0)
    #data = data.apply(lambda row: row.fillna(row.mean()), axis=1)
    #sample_data = data.sample(n=1000, replace=False)
    # for func in [calc_udp_poisson, calc_udp_nbm, calc_udp_gmm, calc_udp_norm, calc_udp_gennorm]:
    #    udp, aic = func(sample_data, aic_test=True)
    #    print(f'Function: {func.__name__}, aic: {aic}')
    #   Function: calc_udp_poisson, aic: 207.66830882318942
    #   Function: calc_udp_nbm, aic: 240.239110524907
    #   Function: calc_udp_gmm, aic: 162.60641335671778
    #   Function: calc_udp_norm, aic: 169.15685593656568
    #   Function: calc_udp_gennorm, aic: 160.8001788887676
    calc_udp_multi_process(False)

# Which one is visually better?
#plt.hist(X, bins=20, normed=True)
#plt.plot(range(260), ss.nbinom.pmf(range(260), np.round(P1[0]), P1[1], np.round(P1[2])), 'g-')
#plt.plot(range(260), ss.nbinom.pmf(range(260), np.round(P2[0]), P2[1], np.round(P2[2])), 'r-')
# Plotting fit vs. RMA values.
#row = udp.iloc[0:,]
#plt.hist(sample_data.iloc[0,:], bins=20, normed=True)
# plt.show()
#plt.hist(row, bins=20, normed=False)
# plt.show()
