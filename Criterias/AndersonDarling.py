from scipy import stats
import numpy as np
import scipy.integrate as integrate
from math import sqrt, pi
from mpmath import nsum, inf, exp, gamma
from IHaveStatistic import IHaveSatistic


class AndersonDarlingCriteria:
    #сравнение с теоретическим: X2 - 'norm'
    def SciPyResult(self, X1, X2):
        return stats.anderson(X1, X2)

    def Result2SamplesBykSamp(self, X1, X2):
        return stats.anderson_ksamp([X1, X2])

    def Result2Samples(self ,X1, X2):
        m = len(X1)
        n = len(X2)
        N = m + n
        X1.sort()
        X2.sort()
        dataAll = np.concatenate((X1, X2))
        dataAll.sort()
        dataAll = dataAll.tolist()

        sum = 0
        sumOfX1 = 0
        for i in range(len(dataAll)-1):
            elem = dataAll[i]
            for j in range(sumOfX1, m):
                if (X1[j] > elem):
                    break
                else:
                    sumOfX1 += 1

            sum += pow((sumOfX1*N - m*(i+1)), 2) / ((i+1)*(N-i-1))

        return IHaveSatistic(sum/(m*n))

    @staticmethod
    def GetStatisticDistribution(stats):
        def a2(jj):
            j = float(jj)
            temp = pow((4.0*j + 1), 2)
            temp2 = 8*stat
            res = gamma(j + 0.5)*(4.0*j + 1)/(gamma(0.5)*gamma(j + 1.0))
            res *= exp(-temp*pi*pi/temp2)
            res *= integrate.quad(lambda y: exp((stat / (8*(y*y+1))) - temp*pi*pi*y*y/temp2), 0, np.inf)[0]
            return pow(-1, j) * res

        result = []
        for stat in stats:
            sum = float(nsum(lambda j: a2(j), [0, inf]))

            st = sqrt(2*pi) / stat
            result.append(st*sum)
        return result
