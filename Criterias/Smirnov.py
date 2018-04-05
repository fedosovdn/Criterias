from scipy import stats

class SmirnovCriteria:
    def SciPyResult(self, X1, X2):
        return stats.kstest(X1, X2)

    def SciPyResult2Samples(self, X, Y):
        return stats.ks_2samp(X, Y);