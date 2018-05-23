from enum import Enum
import scipy.stats as stats
from pandas import DataFrame as df
from pandas import Series as ser
import numpy as np

class CriteriaSide(Enum):
    left = 1
    right = 2
    both = 3

class PowerCalculateHelper:
    @staticmethod
    def CalcualteStats(N, criteria):
        n = 100
        m = 100

        SH0 = []
        SH1 = []

        for i in range(N):
            x1_H0 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H0 = stats.norm.rvs(loc=0, scale=1, size=m)
            SH0.append(criteria.Result2Samples(x1_H0, x2_H0).statistic)

            x1_H1 = stats.norm.rvs(loc=0, scale=1, size=n)
            x2_H1 = stats.norm.rvs(loc=0.1, scale=1, size=m)
            SH1.append(criteria.Result2Samples(x1_H1, x2_H1).statistic)

        print(PowerCalculateHelper.CalculatePower(SH0, SH1, [0.01, 0.05]))

    @staticmethod
    def CalculatePower(statsH0, statsH1,  alphas, criteriaSide = None):
        powers = []
        # quantiles = df(statsH0).quantile(np.ones(len(alphas)) - alphas).values
        quantiles = ser(statsH0).quantile(np.ones(len(alphas)) - alphas).values

        return powers