from scipy import stats

class SmirnovCriteria:
    def SciPyResult(self, X1, X2):
        return stats.kstest(X1, X2)

    def Result2Samples(self, X1, X2):
        return stats.ks_2samp(X1, X2)