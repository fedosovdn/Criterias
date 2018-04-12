import numpy as np
from Helper import Helper
from IHaveStatistic import IHaveSatistic
from scipy.special import iv, gamma
from mpmath import nsum, inf, exp, sqrt

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

        stat = 1/(m*n*(sum0)) * (n * sum2 + m * sum1) - (4*m*n - 1) / (6 * sum0)
        return IHaveSatistic(stat)

    ##Индекс повторяющихся элементов
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
            sum = nsum(lambda j: a1(j), [0, inf])

            st = (1 / sqrt(2 * stat))
            result.append(st*sum)
        return result

    @staticmethod
    def a1(stats):
        v = 0.25

        result = []
        for s in stats:
            a1 = 0
            for i in range(1000):
                z = (4 * i + 1) * (4 * i + 1) / 16.0 / s
                tmp = math.gamma(i + 0.5) / math.gamma(0.5) / math.gamma(i + 1)
                tmp = tmp * math.exp(-z) * math.sqrt(4 * i + 1)
                tmp = tmp * (iv(-v, z) - iv(v, z))
                if (tmp < 0.000000000001):
                    break
                a1 = a1 + tmp
                a1 = a1 / math.sqrt(2 * s)
            result.append(a1)
        return result