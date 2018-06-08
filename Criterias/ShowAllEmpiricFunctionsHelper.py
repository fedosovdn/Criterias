from scipy.stats import norm
from scipy.stats import kstest
from Helper import Helper
from scipy.stats import kstwobign
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
import LehmanRosenblatt as lr
import Smirnov as sm
import AndersonDarling as ad
from statistics import median


class ShowAllEmpiricFunctionsHelper:
    # @staticmethod
    # def ShowAllEmpiricFunctionsOfSmirnovCriteria(n, m, N):
    #     criteria = sm.SmirnovCriteria()
    #     en = np.sqrt(n * m / float(n + m))
    #     func = lambda x: kstwobign.cdf((en + 0.12 + 0.11 / en) * x)
    #     plt = ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'K-S', func)
    #     plt.show()

    @staticmethod
    def ShowAllEmpiricFunctionsOfSmirnovCriteria(n, m, N, loc=0, scale=1):
        criteria = sm.SmirnovCriteria()
        en = np.sqrt(n * m / (n + m))
        func = lambda x: kstwobign.cdf((en + 0.12 + 0.11 / en) * x)
        ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'K-S', func, loc, scale)
        # ShowAllEmpiricFunctionsHelper.GetKolmogorovDistancesOfCriteria(n, m, N, criteria, func, loc, scale)

    @staticmethod
    def ShowAllEmpiricFunctionsOfLehmRosCriteria(n, m, N, loc=0, scale=1):
        criteria = lr.LehmanRosenblattCriteria()
        func = lambda x: lr.LehmanRosenblattCriteria.GetStatisticDistribution(x)
        # ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'L-R', func, loc, scale)
        ShowAllEmpiricFunctionsHelper.GetKolmogorovDistancesOfCriteria(n, m, N, criteria, func, loc, scale)

    @staticmethod
    def ShowAllEmpiricFunctionsOfADCriteria(n, m, N, loc=0, scale=1):
        criteria = ad.AndersonDarlingCriteria()
        func = lambda x: ad.AndersonDarlingCriteria.GetStatisticDistribution(x)
        ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'A-D', func, loc, scale)
        # ShowAllEmpiricFunctionsHelper.GetKolmogorovDistancesOfCriteria(n, m, N, criteria, func, loc, scale)

    @staticmethod
    def GetKolmogorovDistancesOfCriteria(n, m, N, criteria, cdfValues, loc, scale):
        print(f"n: {n}, m: {m}")
        roundingDigitsCounts = [2]  # количество знаков округления значений выборок
        stats = []
        diffCounts = []  # для вывода различных значений в выборках
        for digit in roundingDigitsCounts:
            stats.append([])
            diffCounts.append([])

        nToDisplayUnique = []  # чтобы вывести количество уникальных значений
        count = 10
        for i in range(count):
            nToDisplayUnique.append(N // count * i)

        for i in range(N):
            x1 = norm.rvs(loc=loc, scale=scale, size=n)
            x2 = norm.rvs(loc=loc, scale=scale, size=m)
            for index, digit in enumerate(roundingDigitsCounts):
                x1_rounded = Helper.RoundingArray(x1, digit)
                x2_rounded = Helper.RoundingArray(x2, digit)
                if i in nToDisplayUnique:
                    diffCounts[index].append(len(set(np.concatenate((x1_rounded, x2_rounded)))))
                stats[index].append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)
        for diffCount in diffCounts:
            print(median(diffCount))

        for index, stat in enumerate(stats):
            # print(f"{index} -ый массив статистик")
            print(f"distance by KS test: {kstest(stat, lambda data: cdfValues(data)).statistic}")

    @staticmethod
    def GetPlotsOfCriteria(n, m, N, criteria, shortName ,cdfValues, loc, scale):
        print(f"n: {n}, m: {m}")
        roundingDigitsCounts = [0]  # количество знаков округления значений выборок
        descriptions = (shortName, '0', '1', '2')
        stats = []
        for digit in roundingDigitsCounts:
            stats.append([])

        for i in range(N):
            x1 = norm.rvs(loc=loc, scale=scale, size=n)
            x2 = norm.rvs(loc=loc, scale=scale, size=m)
            for index, digit in enumerate(roundingDigitsCounts):
                x1_rounded = Helper.RoundingArray(x1, digit)
                x2_rounded = Helper.RoundingArray(x2, digit)
                stats[index].append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)

        # lines = []
        allstats = sum(stats, [])
        x = np.linspace(min(allstats), max(allstats), N)
        # lines.append((x, cdfValues(x)))
        plt.plot(x, cdfValues(x))
        for index, stat in enumerate(stats):
            # print(f"{index} -ый массив статистик")
            ecdf = ECDF(stat)
            # lines.append((ecdf.x, ecdf.y))
            plt.plot(ecdf.x, ecdf.y)
        plt.legend(labels=descriptions)
        plt.ylabel("G(Sc|H)")
        plt.xlabel("Sc")

        plt.show()
