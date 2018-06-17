from scipy import stats

from mpmath import nsum, inf, exp, gamma

class SmirnovCriteria:
    def SciPyResult(self, X1, X2):
        return stats.kstest(X1, X2)

    def Result2Samples(self, X1, X2):
        return stats.ks_2samp(X1, X2)

    @staticmethod
    def SmirnovDist(statistics):
        result = []
        for stat in statistics:
            sum = float(nsum(lambda j: ((-1)**j)*exp(-2*j*j*stat*stat), [-inf, inf]))
            result.append(sum)
        return result