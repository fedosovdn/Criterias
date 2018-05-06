from scipy import stats
import numpy as np


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

        return sum/m*n
