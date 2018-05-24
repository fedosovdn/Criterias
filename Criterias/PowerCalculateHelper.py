from enum import Enum
import scipy.stats as stats
from pandas import Series as ser
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from Helper import Helper

class CriteriaSide(Enum):
    left = 1
    right = 2
    both = 3

class PowerCalculateHelper:
    @staticmethod
    def CalcualteStats(N, criteria):
        n = 500
        m = 500

        SH0 = []
        SH1 = []
        SH2 = []
        SH3 = []
        SH4 = []
        SH5 = []

        for i in range(N):
            x1_H0 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H0 = stats.norm.rvs(loc=0, scale=1, size=m)
            SH0.append(criteria.Result2Samples(x1_H0, x2_H0).statistic)

            x1_H1 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H1 = stats.norm.rvs(loc=0.1, scale=1, size=m)
            SH1.append(criteria.Result2Samples(x1_H1, x2_H1).statistic)

            x1_H2 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H2 = stats.norm.rvs(loc=0.5, scale=1, size=m)
            SH2.append(criteria.Result2Samples(x1_H2, x2_H2).statistic)

            x1_H3 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H3 = stats.norm.rvs(loc=0.0, scale=1.1, size=m)
            SH3.append(criteria.Result2Samples(x1_H3, x2_H3).statistic)

            x1_H4 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H4 = stats.norm.rvs(loc=0.0, scale=1.5, size=m)
            SH4.append(criteria.Result2Samples(x1_H4, x2_H4).statistic)

            x1_H5 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H5 = stats.norm.rvs(loc=0.5, scale=1, size=m)
            SH5.append(criteria.Result2Samples(x1_H5, x2_H5).statistic)

        print(PowerCalculateHelper.CalculatePower(SH0, SH1, [0.1, 0.05, 0.025]))

    @staticmethod
    def CalculatePower(statsH0, statsH1,  alphas, criteriaSide = None):
        quantiles = ser(statsH0).quantile(np.ones(len(alphas)) - alphas).values

        ecdf = ECDF(statsH1)

        possibilites = ecdf(quantiles)

        return np.ones(len(alphas)) - possibilites