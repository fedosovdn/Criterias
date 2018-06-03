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
    def ShowAllEmpiricFunctionsOfSmirnovCriteria(n, m, N):
        criteria = sm.SmirnovCriteria()
        en = np.sqrt(n * m / (n + m))
        func = lambda x: kstwobign.cdf((en + 0.12 + 0.11 / en) * x)
        plt = ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'K-S', func)
        plt.show()#пока не закроем окно, следующая итерация не будет

    @staticmethod
    def ShowAllEmpiricFunctionsOfLehmRosCriteria(n, m, N):
        criteria = lr.LehmanRosenblattCriteria()
        func = lambda x: lr.LehmanRosenblattCriteria.GetStatisticDistribution(x)
        plt = ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'L-R', func)
        # plt.show()

    @staticmethod
    def ShowAllEmpiricFunctionsOfADCriteria(n, m, N):
        criteria = ad.AndersonDarlingCriteria()
        func = lambda x: ad.AndersonDarlingCriteria.GetStatisticDistribution(x)
        plt = ShowAllEmpiricFunctionsHelper.GetPlotsOfCriteria(n, m, N, criteria, 'A-D', func)
        # plt.show()

    @staticmethod
    def GetPlotsOfCriteria(n, m, N, criteria, shortName ,cdfValues):
        print(f"n: {n}, m: {m}")
        roundingDigitsCounts = [1]  # количество знаков округления значений выборок
        #descriptions = (shortName, '-', '0', '1', '2')
        descriptions = (shortName, '0', '1', '2')
        #stats = [[]]
        stats = []
        diffCounts = [] # для вывода различных значений в выборках
        for digit in roundingDigitsCounts:
            stats.append([])
            diffCounts.append([])

        nToDisplayUnique = []#чтобы вывести количество уникальных значений
        count = 10
        for i in range(count):
            nToDisplayUnique.append(N//count*i)

        for i in range(N):
            x1 = norm.rvs(loc=0, scale=50, size=n)
            x2 = norm.rvs(loc=0, scale=50, size=m)
            # if i in nToDisplayUnique:
            #     diffCounts += f"{len(set(np.concatenate((x1, x2))))}\t"
            # statistic = criteria.Result2Samples(x1, x2).statistic
            # stats[0].append(statistic)
            for index, digit in enumerate(roundingDigitsCounts):
                x1_rounded = Helper.RoundingArray(x1, digit)
                x2_rounded = Helper.RoundingArray(x2, digit)
                if i in nToDisplayUnique:
                    diffCounts[index].append(len(set(np.concatenate((x1_rounded, x2_rounded)))))
                #stats[index + 1].append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)
                stats[index].append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)
        for diffCount in diffCounts:
            print(median(diffCount))

        lines = []
        # allstats = sum(stats, [])
        # x = np.linspace(min(allstats), max(allstats), N)
        # lines.append(plt.plot(x, cdfValues(x)))
        for index, stat in enumerate(stats):
            # print(f"{index} -ый массив статистик")
            # ecdf = ECDF(stat)
            # lines.append(plt.plot(ecdf.x, ecdf.y))
            print(f"distance by KS test: {kstest(stat, lambda data: cdfValues(data)).statistic}")
        # plt.legend(labels=descriptions)
        # plt.ylabel(shortName)
        # plt.xlabel(f'n: {n}, m: {m}')

        return plt
