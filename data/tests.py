import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture as GMM

"""
# Testing differences between two models. Comparing AIC from two Gamma distributions to two normal distributions.
aic = pd.Series()
#rma_data = pd.read_csv('src/Pathologist/Input/Lung/GSE29013_RMA.txt', delimiter='\t', index_col=0)
rma_data = pd.read_csv('src/Pathologist/Input/ovarian/OV_RMA.txt', delimiter='\t', index_col=0)
rma_data = rma_data.apply(lambda row: row.fillna(row.mean()), axis=1)
#rma_data.replace(0, np.finfo(np.float32).eps, inplace=True)
#aic_data = pd.read_csv('src/Pathologist/Input/Lung/AIC3a.txt', delimiter = '\t', index_col=0)['GSM718769']
aic_data = pd.read_csv('src/Pathologist/Input/ovarian/AIC.txt', delimiter = '\t', index_col=0)['TCGA-01-0628-11A-01R-0363-07']
#AIC is defined as 2k-ln(L), hence we need to minimize it.
#k is the number of estimated parameters in the model.
i = 0
for index, row in rma_data.iterrows():
    if (i%1000 == 0):
        print('#', end="")
    row = row.sort_values().values.reshape(-1, 1)
    gmm = GMM(n_components=2, covariance_type='full', random_state=0).fit(row)
    aic = aic.append(pd.Series([gmm.aic(row)]))
    i += 1

aic.index = aic_data.index
merged_aic = pd.DataFrame({'matlab':aic_data, 'python':aic}, index=aic.index)
merged_aic['delta'] = -aic_data - aic #How much python is smaller that matlab.
merged_aic['ratio'] = merged_aic.delta/(-merged_aic.matlab)
print(merged_aic.ratio.mean()) #-0.028858356977818217  #1.3961894089053348

#Test activity and consistency
output = pd.read_csv('src/Pathologist/Input/Lung/output.csv', index_col=None)
output['pathName'] = output.pathName.str.lower()
matactivity = pd.read_csv('src/Pathologist/Input/Lung/GSE29013_Act.txt', delimiter='\t', encoding='utf-8')
matactivity.columns = ['pathName'] + list(range(0,len(matactivity.columns)-1))
matactivity = pd.melt(matactivity, id_vars=['pathName'], var_name='sampleID') #.sort_values(['Pathway Name', 'sample'])
matactivity.columns = ['pathName', 'sampleID', 'matactivity']
matactivity['sampleID'] = matactivity.sampleID.astype(int)
joined_ac = pd.merge(output, matactivity, how='left')
joined_ac['diff_activity'] = np.abs(joined_ac.Activity - joined_ac.matactivity)
joined_ac.to_csv('src/Pathologist/Input/Lung/joined_ac1.csv', index = None)

################################## Using scipy.stats to fit a distribution.
import matplotlib.pyplot as plt
import scipy
import scipy.stats
import seaborn as sns

size = 30000
x = scipy.arange(size)
y = scipy.int_(scipy.round_(scipy.stats.vonmises.rvs(5,size=size)*47))
n, bins, patches = plt.hist(y, bins=range(48)) #, color='w'
plt.xlabel('Values')
plt.ylabel('Count')
plt.title('Vonmises Distribution')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
#plt.grid(True)

dist_names = ['gamma'] #, 'beta', 'rayleigh', 'norm', 'pareto']

for dist_name in dist_names:
    dist = getattr(scipy.stats, dist_name)
    param = dist.fit(y)
    pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1]) * size
    plt.plot(pdf_fitted, label=dist_name)
    plt.xlim(0,47)
plt.legend(loc='upper right')
plt.show()

################### Pomegranate
import numpy as np
import pandas as pd
from pomegranate import *
from scipy import stats
from sklearn.mixture import GaussianMixture as GMM

X=np.array([2.44515009, 5.50218787, 2.78317919, 3.57474729, 2.28051024,
        2.30735988, 3.14827695, 4.66004906, 2.17413335, 1.03239638,
        4.49822428, 4.02841676, 1.57646796, 2.9348154, 3.45762544,
        1.86445921, 2.05579539, 1.01855821, 3.18763872, 1.58823228,
        2.5070555, 2.33347747, 1.79635148, 4.18190807, 3.79100177,
        2.0217273, 3.56959567, 1.74117038, 3.41538012, 5.36599237,
        4.50731341, 2.28592422, 3.69336315, 2.47597859, 0.93218948,
        1.91021177, 3.78694875, 2.30966696, 3.30454069, 3.73907217,
        2.29007995, 1.49829356, 2.8240885, 2.22911464, 2.8449105,
        2.83521877, 2.57443547, 1.86617824, 2.98622416, 1.9268491,
        2.33751329, 2.96892834, 4.32908499, 2.5824886, 2.29107388])

gmm = GMM(n_components=2, covariance_type='full', random_state=0).fit(X.reshape(-1,1))
labels = gmm.predict(X.reshape(-1,1))
data = pd.read_csv('src/Pathologist/Input/Lung/GSE29013_RMA.txt', delimiter='\t', index_col=0)
data.replace(0, np.finfo(np.float32).eps, inplace=True)
d1 = ExponentialDistribution(5, 2)
d2 = ExponentialDistribution(4, 1)
#d1 = distributions.IndependentComponentsDistributions([distributions.NormalDistribution(5, 2), ExponentialDistribution(1), LogNormalDistribution(0.4, 0.1)])
#d1 = IndependentComponentsDistribution([GammaDistribution(3, 1), GammaDistribution(2, 3)])
model = GeneralMixtureModel([d1, d2], weights=[0.66, 0.34])
X=data.T[['AFFX-TrpnX-M_at']].values.reshape(1,-1)
model = GeneralMixtureModel.from_samples(NormalDistribution, n_components=2, X=X)
#model = NormalDistribution.from_samples(X)


################ https://www.kaggle.com/charel/learn-by-example-expectation-maximization
class GeneralMixture:
    "Model mixture of two univariate Gaussians and their EM estimation"
    def __init__(self, data, mix=.5):
        self.data = data
        # todo the Algorithm would be numerical enhanced by normalizing the data first, next do all the EM steps and do the de-normalising at the end
        # init with multiple gaussians
        self.one = stats.gamma(1, 10)
        self.two = stats.gamma(1, 9)
        # as well as how much to mix them
        self.mix = mix

    def Estep(self):
        "Perform an E(stimation)-step, assign each point to gaussian 1 or 2 with a percentage"
        # compute weights
        self.loglike = 0.  # = log(p = 1)
        for datum in self.data:
            # unnormalized weights
            wp1 = self.one.pdf(datum) * self.mix
            wp2 = self.two.pdf(datum) * (1. - self.mix)
            # compute denominator
            den = wp1 + wp2
            # normalize (numerator/denominator)
            wp1 /= den
            wp2 /= den  # wp1+wp2= 1, it either belongs to distribution 1 or distribution 2
            # add into loglike
            self.loglike += np.log(wp1 + wp2)  # freshening up self.loglike in the process
            # yield weight tuple
            yield (wp1, wp2)

    def Mstep(self, weights):
        "Perform an M(aximization)-step"
        # compute denominators
        (left, rigt) = zip(*weights)
        one_den = sum(left)
        two_den = sum(rigt)

        # compute new means
        self.one.mu = sum(w * d for (w, d) in zip(left, self.data)) / one_den
        self.two.mu = sum(w * d for (w, d) in zip(rigt, self.data)) / two_den

        # compute new sigmas
        self.one.sigma = np.sqrt(sum(w * ((d - self.one.mu) ** 2)
                                  for (w, d) in zip(left, self.data)) / one_den)
        self.two.sigma = np.sqrt(sum(w * ((d - self.two.mu) ** 2)
                                  for (w, d) in zip(rigt, self.data)) / two_den)
        # compute new mix
        self.mix = one_den / len(self.data)

    def iterate(self, N=1, verbose=False):
        "Perform N iterations, then compute log-likelihood"
        for i in range(1, N + 1):
            self.Mstep(self.Estep())  # The heart of the algorithm, perform E-step and next M-step
            if verbose:
                print('{0:2} {1}'.format(i, self))
        self.Estep()  # to freshen up self.loglike

    def pdf(self, x):
        return (self.mix) * self.one.pdf(x) + (1 - self.mix) * self.two.pdf(x)

    def __repr__(self):
        return (f'GeneralMixture({self.one}, {self.two}, mix={self.mix}')

    def __str__(self):
        return (f'Mixture: {self.one}, {self.two}, mix={self.mix}')


# See the algorithem in action
n_iterations = 5
best_mix = None
best_loglike = float('-inf')
gm = GeneralMixture(X)
for _ in range(n_iterations):
    try:
        # train!
        gm.iterate(verbose=True)
        if gm.loglike > best_loglike:
            best_loglike = gm.loglike
            best_mix = gm.mix

    except (ZeroDivisionError, ValueError, RuntimeWarning):  # Catch division errors from bad starts, and just throw them out...
        print("Exception. one less")
        pass
print(f"Done. Best mixture:{best_mix}")

################ https://gist.github.com/rmcgibbo/7318011
import numpy as np
import pandas as pd
import scipy.special
from scipy.misc import logsumexp
import matplotlib.pyplot as pp
import scipy.stats

def invpsi(y):
    #Inverse digamma (psi) function.  The digamma function is the
    #derivative of the log gamma function.
    # Adapted from matlab code in PMTK (https://code.google.com/p/pmtk3), copyright
    # Kevin Murphy and available under the MIT license.

    # Newton iteration to solve digamma(x)-y = 0
    x = np.exp(y) + 0.5
    mask = y < -2.22
    x[mask] = 1.0 / (y[mask] - scipy.special.psi(1))

    # never more than 5 iterations required
    for i in range(5):
        x = x - (scipy.special.psi(x) - y) / scipy.special.polygamma(1, x)
    return x


def fit(x, k=3):
    alpha = 10 * np.random.rand(k)
    rate = 10 * np.random.rand(k)
    pi = np.ones(k) #Weights of distributions.
    log_x = np.log(x)

    for i in range(50):
        # log probability of each data point in each component
        #gammaln - log of the absolute value of gamma function for real x.
        logg = alpha * np.log(rate) - scipy.special.gammaln(alpha) + \
               np.multiply.outer(log_x, alpha - 1) - np.multiply.outer(x, rate)
        logp = np.log(pi) + logg - logsumexp(logg, axis=1, b=pi[np.newaxis])[:, np.newaxis]
        p = np.exp(logp)

        # new mixing weights
        pi = np.mean(np.exp(logp), axis=0)

        # new rate and scale parameters
        A = np.einsum('i,ij->j', log_x, p)
        B = np.einsum('j,ij->j', np.log(rate), p)
        alpha_argument = (A + B) / np.sum(p, axis=0)
        rate = alpha * np.sum(p, axis=0) / np.einsum('i,ij->j', x, p)
        # when the fit is bad (early iterations), this conditional maximum
        # likelihood update step is not guarenteed to keep alpha positive,
        # which causes the next iteration to be f*cked.
        alpha = np.maximum(invpsi(alpha_argument), 1e-8)

    x = np.linspace(0.001, np.max(x), 1000)
    g = (np.power(rate, alpha) / scipy.special.gamma(alpha)) * np.power.outer(x, alpha - 1) * \
         np.exp(-np.multiply.outer(x, rate))
    ax2 = pp.gca().twinx()
    ax2.plot(x, g[:, 0])
    ax2.plot(x, g[:, 1])

x = np.concatenate((scipy.stats.distributions.gamma(9, 3).rvs(500),
                    scipy.stats.distributions.gamma(3, 0.2).rvs(500)))
pp.hist(x, bins=50, alpha=0.3)
fit(x, k=2)
pp.show()

"""
"""
cmplx = cmplx.loc[cmplx[1] != 'rna'] # RNA complexes never have links.
#cmplx_1 = cmplx.loc[cmplx[1] != 'complex'].copy()
#In case of protein its ID doesn't matter only the links.
#cmplx_1.loc[cmplx_1[1] == 'protein', 'pr'] = cmplx[4].apply(lambda x: max([parse_link(i) for i in str(x).split(',')]))
#cmplx_1.loc[cmplx_1[1] == 'compound', 'pr'] = 1

#cmplx_basic = cmplx_basic.rename(columns={0:3}) #Columns are now (3==complex ID, 'pr').
##cmplx_non_basic = cmplx.loc[cmplx[1] == 'complex'].copy()
##cmplx_non_basic = pd.merge(cmplx_non_basic, cmplx_basic, on=3, how='left') #By merging we can find non-basic complexes UDP.
#How to handle cases where basic complexes are not found?
##cmplx_non_basic['pr'] = cmplx_non_basic['pr'].fillna(1)
#Calculate the UDP for a complex molecule now that we have UDP for all the ingredients.
##cmplx_non_basic = cmplx_non_basic[[0, 'pr']].groupby(0, as_index=False).prod()
##cmplx_non_basic = cmplx_non_basic.rename(columns={0:3})
"""
"""
################## Plotting mixture of normal distributions.
for index, row in rma_data.iterrows():
    row = row.sort_values().values.reshape(-1, 1)
    if (index == '1552327_at'):
        print('1552327_at')
        break
    gmm = GMM(n_components=2, covariance_type='full', random_state=0).fit(row)
    aic = aic.append(pd.Series([gmm.aic(row)]))
#n, bins, patches = plt.hist(row.values, bins=range(22))  # , color='w'
#pred = gmm.predict_proba(row.values.reshape(-1, 1))[:,1]
#plt.plot(pred, label='GaussianMixture')
#plt.xlim(0, 47)
#plt.legend(loc='upper right')
#labels = gmm.predict(row.values.reshape(-1,1))
#https://stackoverflow.com/questions/24878729/how-to-construct-and-plot-uni-variate-gaussian-mixture-using-its-parameters-in-p
x = np.linspace(min(row), max(row), 50)
means = gmm.means_
stdevs = np.sqrt(gmm.covariances_)[:,0]
weights = gmm.weights_
pdfs = [p * stats.norm.pdf(x, mu, sd) for mu, sd, p in zip(means, stdevs, weights)]
density = np.sum(np.array(pdfs), axis=0)
plt.plot(x, density)
plt.show()

########## Plotting mixture of Gamma distributions.
#alphas = [591.1995546, 331.5439896]
#betas = [0.093856608, 0.02209614]
#weights = [0.83, 0.17]
first = [56.6074694, 0.01417483]
second = [8.10133829, 0.424369654]
weights = [0.596102493, 0.403897507]
x = np.linspace(0, 9, 100)
alphas = [56.6074694, 8.10133829]
betas = [0.01417483, 0.424369654]
pdfs = [p * stats.gamma.pdf(x, a=mu, scale=sd) for mu, sd, p in zip(alphas, betas, weights)]
density = np.sum(np.array(pdfs), axis=0)
#density = 0.59*stats.gamma.pdf(x, a=first[0], scale=first[1]) + 0.41*stats.gamma.pdf(x, a=second[0], scale=second[1])
plt.plot(x, density) # , weights[0]*first[0]+weights[0]*second[0]
plt.show()

#https://www.datacamp.com/community/tutorials/probability-distributions-python
#This does not take into account the weights.
data_0 = stats.gamma.rvs(a=first[0], scale=first[1], size=1000)
data_1 = stats.gamma.rvs(a=second[0], scale=second[1], size=1000)
sns.distplot(data_0,kde=True,bins=100,color='skyblue',hist_kws={"linewidth": 15,'alpha':1})
sns.distplot(data_1,kde=True,bins=100,color='blue',hist_kws={"linewidth": 15,'alpha':1})
plt.xlabel('Gamma Distribution')
plt.ylabel('Frequency')
plt.show()

######
red = 255; #i.e. FF
green = 0;
stepSize = 4 #how many colors do you want?
while(green < 255):
    green += stepSize
    if(green > 255):
        green = 255
    print(str(hex(red)), str(hex(green)), 0) #assume output is function that takes RGB

while(red > 0):
    red -= stepSize;
    if(red < 0):
        red = 0
    print(str(hex(red)), str(hex(green)), 0) #assume output is function that takes RGB

red = 0xFF0000
colors = ['00FF00','11FF00','22FF00','33FF00','44FF00','55FF00','66FF00','77FF00','88FF00',
          '99FF00','AAFF00','BBFF00','CCFF00','DDFF00','EEFF00','FFFF00','FFEE00','FFDD00',
          'FFCC00','FFBB00','FFAA00','FF9900','FF8800','FF7700','FF6600','FF5500','FF4400',
          'FF3300','FF2200','FF1100','FF0000']
######
"""
#data = pd.read_csv('c:/users/Admin/Downloads/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', delimiter='\t', index_col = 0, error_bad_lines=False)
data = pd.read_csv('c:/users/Admin/Downloads/ESCA.rnaseqv2.txt', delimiter='\t', index_col = 0, error_bad_lines=False)
data1 = data[data.columns[::3]]
data1.to_csv('c:/users/Admin/Downloads/ESCA.csv')

###### Remove Entre gene ID.
data = pd.read_csv('c:/users/Admin/Downloads/ESCA.csv', index_col=0)
#data['index'] = data['index'].apply(lambda x: x.split('|')[0])
data.reset_index(inplace=True)
data['index'] = data['index'].apply(lambda x: x.split('|')[0])
data = data.loc[data['index'] != '?']
data.index = data['index']
data.drop('index', axis=1, inplace=True)
data.to_csv('c:/users/Admin/downloads/ESCA_clean.csv')

######
#Use statsmodels for fitting NegativeBinomial.
import numpy as np
import pandas as pd
from statsmodels.discrete.discrete_model import NegativeBinomialP
from math import log

data = pd.read_csv('ACC.rnaseq.csv', index_col=0)
data = data.apply(lambda row: row.fillna(row.mean()), axis=1)
exog = np.ones(len(data.columns))
my_udp = np.empty((0, len(data.columns)))
row = data.loc['ABCA6|23460']
row = row.astype(int).values #.reshape(-1,1)
if(max(row)<5 or np.count_nonzero(row)<10): #(row == 0).all() or np.sum(row)<200
   my_udp = np.append(my_udp, [np.zeros(len(data.columns))], axis=0)
row[np.argsort(row)[-3]] = row[np.argsort(row)[-2]] = row[np.argsort(row)[-1]] = row[np.argsort(row)[-4]]
model_nb = NegativeBinomialP(row, exog, p=2)
res_nb = model_nb.fit(method='bfgs', maxiter=5000, maxfun=5000, disp=0)

# Copyright (C) 2014 Gokcen Eraslan
import warnings
import sys

with warnings.catch_warnings():
   np.seterr(all='raise')
   try:
      result = np.sum(gammaln(X + r)) \
               - np.sum(np.log(factorial(X))) \
               - N * (gammaln(r)) \
               + N * r * np.log((1 - p if p < 1 else infinitesimal)) \
               + np.sum(X * np.log(p))
      # Both Wikipedia and np.negative_binomial have k and r swapped.
      # + N*r*np.log(p) \
      # + np.sum(X*np.log(1-(p if p < 1 else 1-infinitesimal)))
   except Warning:
      print(p)
      sys.stdout.flush()
      sys.exit()

######
#Colors from matplotlib
cnames = {'aliceblue':'#F0F8FF','antiquewhite':'#FAEBD7','aqua':'#00FFFF','aquamarine':'#7FFFD4','azure':'#F0FFFF','beige':'#F5F5DC',
'bisque':'#FFE4C4','black':'#000000','blanchedalmond':'#FFEBCD','blue':'#0000FF','blueviolet':'#8A2BE2','brown':'#A52A2A','burlywood':'#DEB887',
'cadetblue':'#5F9EA0','chartreuse':'#7FFF00','chocolate':'#D2691E','coral':'#FF7F50','cornflowerblue':'#6495ED','cornsilk':'#FFF8DC','crimson':'#DC143C',
'cyan':'#00FFFF','darkblue':'#00008B','darkcyan':'#008B8B','darkgoldenrod':'#B8860B','darkgray':'#A9A9A9','darkgreen':'#006400',
'darkkhaki':'#BDB76B','darkmagenta':'#8B008B','darkolivegreen':'#556B2F','darkorange':'#FF8C00','darkorchid':'#9932CC','darkred':'#8B0000',
'darksalmon':'#E9967A','darkseagreen':'#8FBC8F','darkslateblue':'#483D8B','darkslategray':'#2F4F4F','darkturquoise':'#00CED1','darkviolet':'#9400D3',
'deeppink':'#FF1493','deepskyblue':'#00BFFF','dimgray':'#696969','dodgerblue':'#1E90FF','firebrick':'#B22222','floralwhite':'#FFFAF0',
'forestgreen':'#228B22','fuchsia':'#FF00FF','gainsboro':'#DCDCDC','ghostwhite':'#F8F8FF','gold':'#FFD700','goldenrod':'#DAA520',
'gray':'#808080','green':'#008000','greenyellow':'#ADFF2F','honeydew':'#F0FFF0','hotpink':'#FF69B4','indianred':'#CD5C5C','indigo':'#4B0082',
'ivory':'#FFFFF0','khaki':'#F0E68C','lavender':'#E6E6FA','lavenderblush':'#FFF0F5','lawngreen':'#7CFC00','lemonchiffon':'#FFFACD',
'lightblue':'#ADD8E6','lightcoral':'#F08080','lightcyan':'#E0FFFF','lightgoldenrodyellow':'#FAFAD2','lightgreen':'#90EE90','lightgray':'#D3D3D3',
'lightpink':'#FFB6C1','lightsalmon':'#FFA07A','lightseagreen':'#20B2AA','lightskyblue':'#87CEFA','lightslategray':'#778899','lightsteelblue':'#B0C4DE',
'lightyellow':'#FFFFE0','lime':'#00FF00','limegreen':'#32CD32','linen':'#FAF0E6','magenta':'#FF00FF','maroon':'#800000','mediumaquamarine':'#66CDAA',
'mediumblue':'#0000CD','mediumorchid':'#BA55D3','mediumpurple':'#9370DB','mediumseagreen':'#3CB371','mediumslateblue':'#7B68EE',
'mediumspringgreen':'#00FA9A','mediumturquoise':'#48D1CC','mediumvioletred':'#C71585','midnightblue':'#191970','mintcream':'#F5FFFA',
'mistyrose':'#FFE4E1','moccasin':'#FFE4B5','navajowhite':'#FFDEAD','navy':'#000080','oldlace':'#FDF5E6','olive':'#808000',
'olivedrab':'#6B8E23','orange':'#FFA500','orangered':'#FF4500','orchid':'#DA70D6','palegoldenrod':'#EEE8AA','palegreen':'#98FB98',
'paleturquoise':'#AFEEEE','palevioletred':'#DB7093','papayawhip':'#FFEFD5','peachpuff':'#FFDAB9','peru':'#CD853F','pink':'#FFC0CB',
'plum':'#DDA0DD','powderblue':'#B0E0E6','purple':'#800080','red':'#FF0000','rosybrown':'#BC8F8F','royalblue':'#4169E1',
'saddlebrown':'#8B4513','salmon':'#FA8072','sandybrown':'#FAA460','seagreen':'#2E8B57','seashell':'#FFF5EE','sienna':'#A0522D',
'silver':'#C0C0C0','skyblue':'#87CEEB','slateblue':'#6A5ACD','slategray':'#708090','snow':'#FFFAFA','springgreen':'#00FF7F',
'steelblue':'#4682B4','tan':'#D2B48C','teal':'#008080','thistle':'#D8BFD8','tomato':'#FF6347','turquoise':'#40E0D0','violet':'#EE82EE',
'wheat':'#F5DEB3','white':'#FFFFFF','whitesmoke':'#F5F5F5','yellow':'#FFFF00','yellowgreen':'#9ACD32'}

# udp = pd.read_csv(udp_file, index_col=0)
# udp = udp.iloc[:,sample_num]
# series = [udp, udp]
# df = pd.concat(series, axis=1)
# df.to_csv(f'data/output_udp_sample_{sampleID}.csv')
# We calculate the UDPs twice since calc_cativity_consistency accepts only datframes.

<div class="content-section">
<h3>Logs</h3>
<textarea rows="4" cols="30">
    {% with messages = get_flashed_messages(category_filter=["logs"]) %}
    {% if messages %}
            <ul class=flashes>
            {% for message in messages %}
                <li>{{ message }}</li>
            {% endfor %}
            </ul>
    {% endif %}
    {% endwith %}
</textarea>
</div>

######
# Filter only probes that appear in paths. We collect all links that are mentioned in either the paths or complexes database.
def filter_probes(self):
    path_links = set()

    # Find all links that are used in paths' proteins.
    self.orig_paths.molLink.loc[self.orig_paths.molType == 'protein'].apply(
        lambda x: [path_links.add(i) for i in str(x).split(',')])

    # Find all links that are used in paths' complexes.
    path_cmplx = self.orig_paths.loc[self.orig_paths.molType == 'complex'].copy()
    cmplx_basic = self.orig_cmplx[['complex_ID', 'molLink']].copy().dropna(
        subset=['molLink'])  # non-basic complexes don't have links.
    cmplx_basic.columns = ['molID', 'molLink']
    cmplx_non_basic = self.orig_cmplx[['complex_ID', 'molID', 'molLink']].copy()
    cmplx_non_basic = cmplx_non_basic[cmplx_non_basic['molLink'].isnull()]  # non-basic complexes don't have links.
    cmplx_non_basic.drop(['molLink'], inplace=True, axis=1)
    # Add links to the non-basic complexes.
    cmplx_non_basic = pd.merge(cmplx_non_basic, cmplx_basic, how='left')  # complex_ID, molID, molLink
    cmplx_non_basic.drop(['molID'], inplace=True, axis=1)  # No need for the basic complex ID.
    # Merge the two type of complexes.
    cmplx_non_basic.columns = ['molID', 'molLink']
    cmplx_list = pd.concat([cmplx_basic, cmplx_non_basic], axis=0).reset_index(drop=True)
    cmplx_list.columns = ['molNum', 'molLink']

    # Merge with paths.
    path_cmplx.drop(['molLink'], inplace=True, axis=1)
    path_cmplx = pd.merge(path_cmplx, cmplx_list, how='left')
    path_cmplx.apply(lambda x: [path_links.add(i) for i in str(x).split(',')])

    new_index = pd.Series()
    print(path_links)