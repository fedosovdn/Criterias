import numpy as np
from Helper import Helper
from IHaveStatistic import IHaveSatistic
from scipy.special import iv
from mpmath import nsum, inf, exp, gamma
from math import sqrt

import math


class LehmanRosenblattCriteria:
    def Result2Samples(self, X1, X2):
        m = len(X1)
        n = len(X2)
        sum0 = m + n
        X1.sort()
        X2.sort()
        dataAll = np.concatenate((X1, X2))
        dataAll.sort()
        dataAll = dataAll.tolist()

        # X1 = X1.tolist()
        # X2 = X2.tolist()
        sum1 = 0
        for i in range(m):
            elem = X1[i]
            sum1 += pow(self.rangOfRepeatingElems(dataAll, elem) - i, 2)

        sum2 = 0
        for i in range(n):
            elem = X2[i]
            sum2 += pow(self.rangOfRepeatingElems(dataAll, elem) - i, 2)

        stat = (n * sum2 + m * sum1)/(m*n*sum0) - (4*m*n - 1) / (6 * sum0)
        return IHaveSatistic(stat if stat != 0 else 1E-15)#иногда статистика получается равна 0, не должно быть такого
        #return IHaveSatistic(stat)

    #Индекс повторяющихся элементов
    def rangOfRepeatingElems(self, variationalSeries, value):
        first = variationalSeries.index(value)
        last = Helper.GetLastIndexOf(variationalSeries, value)
        return (first + last) / 2

    @staticmethod
    def GetStatisticDistribution(stats):
        def a1(jj):
            j = float(jj)
            el = (4 * j + 1)*(4 * j + 1) / 16.0 / stat
            temp = gamma(j + 0.5)*sqrt(4.0*j + 1)/(gamma(0.5)*gamma(j + 1.0))
            bessel = (iv(-0.25, el) - iv(0.25, el))
            return temp * exp(-el) * bessel if not math.isnan(bessel) else 0.0

        result = []
        for stat in stats:
            sum = float(nsum(lambda j: a1(j), [0, inf]))

            st = (1 / sqrt(2 * stat))
            result.append(st*sum)
        return result
