import numpy as np
from Helper import Helper
from IHaveStatistic import IHaveSatistic


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

    def rangOfRepeatingElems(self, variationalSeries, value):
        first = variationalSeries.index(value)
        last = Helper.GetLastIndexOf(variationalSeries, value)
        return (first + last) / 2