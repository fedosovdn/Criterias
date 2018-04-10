from scipy import stats


class AndersonDarlingCriteria:
    #сравнение с теоретическим: X2 - 'norm'
    def SciPyResult(self, X1, X2):
        return stats.anderson(X1, X2)

    def Result2Samples(self, X1, X2):
        return stats.anderson_ksamp([X1, X2])