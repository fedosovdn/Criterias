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
    def CalcualteStats(n, m, N, criteria, digit):
        SH0 = []
        SH1 = []
        SH2 = []
        SH3 = []
        SH4 = []
        SH5 = []

        for i in range(N):
            rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            x1_H0 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0, scale=1, size=m)
            x2_H0 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH0.append(criteria.Result2Samples(x1_H0, x2_H0).statistic)

            rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            x1_H1 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0.1, scale=1, size=m)
            x2_H1 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH1.append(criteria.Result2Samples(x1_H1, x2_H1).statistic)

            rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            x1_H2 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0.5, scale=1, size=m)
            x2_H2 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH2.append(criteria.Result2Samples(x1_H2, x2_H2).statistic)

            rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            x1_H3 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0, scale=1.1, size=m)
            x2_H3 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH3.append(criteria.Result2Samples(x1_H3, x2_H3).statistic)

            rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            x1_H4 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0, scale=1.5, size=m)
            x2_H4 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH4.append(criteria.Result2Samples(x1_H4, x2_H4).statistic)

            rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            x1_H5 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.logistic.rvs(loc=0, scale=1, size=m)
            x2_H5 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH5.append(criteria.Result2Samples(x1_H5, x2_H5).statistic)

        print(PowerCalculateHelper.CalculatePower(SH0, SH1, [0.1, 0.05, 0.025]))
        print(PowerCalculateHelper.CalculatePower(SH0, SH2, [0.1, 0.05, 0.025]))
        print(PowerCalculateHelper.CalculatePower(SH0, SH3, [0.1, 0.05, 0.025]))
        print(PowerCalculateHelper.CalculatePower(SH0, SH4, [0.1, 0.05, 0.025]))
        print(PowerCalculateHelper.CalculatePower(SH0, SH5, [0.1, 0.05, 0.025]))

    @staticmethod
    def CalculatePower(statsH0, statsH1,  alphas, criteriaSide = None):
        quantiles = ser(statsH0).quantile(np.ones(len(alphas)) - alphas).values

        ecdf = ECDF(statsH1)

        possibilites = ecdf(quantiles)

        print()

        return np.ones(len(alphas)) - possibilites