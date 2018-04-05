from scipy import stats


class AndersonDarlingCriteria:
    def SciPyResult(self, X1, X2):
        return stats.anderson(X1, X2)

    def SciPiResultKSamples(self, samples, midrank=True):
        return stats.anderson_ksamp(samples, midrank)