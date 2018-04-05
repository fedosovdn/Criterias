import numpy as np
from Helper import Helper


class LehmanRosenblattCriteria:
    def Result(self, X1, X2):
        m = len(X1)
        n = len(X2)
        sum0 = m + n
        unitedVariationalSeries = np.concatenate((X1, X2))
        unitedVariationalSeries.sort()
        unitedVariationalSeries = unitedVariationalSeries.tolist()
        # X1 = X1.tolist()
        # X2 = X2.tolist()
        sum1 = 0
        for i in range(m):
            elem = X1[i]
            sum1 += pow(self.rangOfRepeatingElems(unitedVariationalSeries, elem) - i, 2)

        sum2 = 0
        for i in range(n):
            elem = X2[i]
            sum2 += pow(self.rangOfRepeatingElems(unitedVariationalSeries, elem) - i, 2)

        stat = 1/(m*n*(sum0)) * (n * sum2 + m * sum1) - (4*m*n - 1) / (6 * sum0)
        return stat

    def rangOfRepeatingElems(self, variationalSeries, value):
        first = variationalSeries.index(value)
        last = Helper.GetLastIndexOf(variationalSeries, value)
        return (first + last) / 2