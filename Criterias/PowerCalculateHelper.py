import scipy.stats as stats
from pandas import Series as ser
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from Helper import Helper
import matplotlib.pyplot as plt
import LehmanRosenblatt as lr

class PowerCalculateHelper:
    @staticmethod
    def CalculateStats(n, m, N, criteria, digit):
        SH0 = []
        SH1 = []
        SH2 = []
        SH3 = []
        SH4 = []
        SH5 = []

        for i in range(N):
            # rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            # x1_H0 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # rvs = stats.norm.rvs(loc=0, scale=1, size=m)
            # x2_H0 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # SH0.append(criteria.Result2Samples(x1_H0, x2_H0).statistic)
            #
            # rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            # x1_H1 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # rvs = stats.norm.rvs(loc=0.1, scale=1, size=m)
            # x2_H1 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # SH1.append(criteria.Result2Samples(x1_H1, x2_H1).statistic)
            #
            # rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            # x1_H2 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # rvs = stats.norm.rvs(loc=0.5, scale=1, size=m)
            # x2_H2 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # SH2.append(criteria.Result2Samples(x1_H2, x2_H2).statistic)
            #
            # rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            # x1_H3 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # rvs = stats.norm.rvs(loc=0, scale=1.1, size=m)
            # x2_H3 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # SH3.append(criteria.Result2Samples(x1_H3, x2_H3).statistic)
            #
            # rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            # x1_H4 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # rvs = stats.norm.rvs(loc=0, scale=1.5, size=m)
            # x2_H4 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # SH4.append(criteria.Result2Samples(x1_H4, x2_H4).statistic)
            #
            # rvs = stats.norm.rvs(loc=0, scale=1, size=n)
            # x1_H5 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # rvs = stats.logistic.rvs(loc=0, scale=1, size=m)
            # x2_H5 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            # SH5.append(criteria.Result2Samples(x1_H5, x2_H5).statistic)



            rvs = stats.norm.rvs(loc=0, scale=2, size=n)
            x1_H4 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0, scale=2, size=m)
            x2_H4 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH4.append(criteria.Result2Samples(x1_H4, x2_H4).statistic)

            rvs = stats.norm.rvs(loc=0, scale=3, size=n)
            x1_H5 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            rvs = stats.norm.rvs(loc=0, scale=3, size=m)
            x2_H5 = rvs if digit == '-' else Helper.RoundingArray(rvs, digit)
            SH5.append(criteria.Result2Samples(x1_H5, x2_H5).statistic)

        # Вычисление мощностей
        # print(PowerCalculateHelper.CalculatePower(SH0, SH1, [0.1, 0.05, 0.025]))
        # print(PowerCalculateHelper.CalculatePower(SH0, SH2, [0.1, 0.05, 0.025]))
        # print(PowerCalculateHelper.CalculatePower(SH0, SH3, [0.1, 0.05, 0.025]))
        # print(PowerCalculateHelper.CalculatePower(SH0, SH4, [0.1, 0.05, 0.025]))
        # print(PowerCalculateHelper.CalculatePower(SH0, SH5, [0.1, 0.05, 0.025]))

        # Отрисовка распределений статситик при альтернативах
        # statistics = [SH0, SH1, SH2, SH3, SH4, SH5]
        # descriptions = ['A-D', 'G(Sc|H0)', 'G(Sc|H1)', 'G(Sc|H2)', 'G(Sc|H3)', 'G(Sc|H4)', 'G(Sc|H5)']
        # func = lambda x: ad.AndersonDarlingCriteria.GetStatisticDistribution(x)

        statistics = [SH4, SH5]
        descriptions = ['L-R', 'G(Sc|scale=2, n=m=200)', 'G(Sc|scale=3, n=m=200)']
        func = lambda x: lr.LehmanRosenblattCriteria.GetStatisticDistribution(x)
        PowerCalculateHelper.ShowPlotsByStats(statistics, descriptions, func, N)

    @staticmethod
    def CalculatePower(statsH0, statsH1,  alphas, criteriaSide = None):
        quantiles = ser(statsH0).quantile(np.ones(len(alphas)) - alphas).values

        ecdf = ECDF(statsH1)

        possibilites = ecdf(quantiles)

        print()

        return np.ones(len(alphas)) - possibilites

    @staticmethod
    def ShowPlotsByStats(statisticsArray, descriptions, cdfValues, N):
        lines = []
        allstats = sum(statisticsArray, [])
        x = np.linspace(min(allstats), max(allstats), N)
        # x = np.linspace(0.1, 30, N)
        lines.append(plt.plot(x, cdfValues(x)))
        for index, stat in enumerate(statisticsArray):
            ecdf = ECDF(stat)
            lines.append(plt.plot(ecdf.x, ecdf.y))
        plt.legend(labels=descriptions)
        plt.ylabel("G(Sc|H)")
        plt.xlabel("Sc")
        plt.show()

    @staticmethod
    def ShowTheorDistrs():
        descriptions = ['F(x|H0)', 'F(x|H1)', 'F(x|H2)', 'F(x|H3)', 'F(x|H4)', 'F(x|H5)']
        # descriptions = ['F(x|H0)',  'F(x|H5)']
        x = np.linspace(-5, 5, 500)
        plt.plot(x, stats.norm.cdf(x, loc=0, scale=1))

        x = np.linspace(-5, 5, 500)
        plt.plot(x, stats.norm.cdf(x, loc=0.1, scale=1))

        x = np.linspace(-5, 5, 500)
        plt.plot(x, stats.norm.cdf(x, loc=0.5, scale=1))

        x = np.linspace(-5, 5, 500)
        plt.plot(x, stats.norm.cdf(x, loc=0, scale=1.1))

        x = np.linspace(-5, 5, 500)
        plt.plot(x, stats.norm.cdf(x, loc=0, scale=1.5))

        x = np.linspace(-5, 5, 500)
        plt.plot(x, stats.logistic.cdf(x, loc=0, scale=1))

        plt.legend(labels=descriptions)
        plt.annotate('frffr', xy=(1, 0.5), xytext=(-3, 1))
        plt.ylabel("F(x|H)")
        plt.xlabel("x")
        plt.show()

    @staticmethod
    def ShowDistrPlots(n, digit):
        samples = []
        rvs = stats.norm.rvs(loc=0, scale=1, size=n)
        samples.append(rvs if digit == '-' else Helper.RoundingArray(rvs, digit))
        # rvs = stats.norm.rvs(loc=0.1, scale=1, size=n)
        # samples.append(rvs if digit == '-' else Helper.RoundingArray(rvs, digit))
        # rvs = stats.norm.rvs(loc=0.5, scale=1, size=n)
        # samples.append(rvs if digit == '-' else Helper.RoundingArray(rvs, digit))
        # rvs = stats.norm.rvs(loc=0.0, scale=1.1, size=n)
        # samples.append(rvs if digit == '-' else Helper.RoundingArray(rvs, digit))
        # rvs = stats.norm.rvs(loc=0.0, scale=1.5, size=n)
        # samples.append(rvs if digit == '-' else Helper.RoundingArray(rvs, digit))
        rvs = stats.logistic.rvs(loc=0.0, scale=1.0, size=n)
        samples.append(rvs if digit == '-' else Helper.RoundingArray(rvs, digit))
        for sample in samples:
            ecdf = ECDF(sample)
            plt.plot(ecdf.x, ecdf.y)
        descriptions = ['F(x|H0)', 'F(x|H5)']
        # descriptions = ['F(x|H0)', 'F(x|H1)', 'F(x|H2)', 'F(x|H3)', 'F(x|H4)', 'F(x|H5)']
        plt.legend(labels=descriptions)
        plt.ylabel("F(x|H)")
        plt.xlabel("x")
        plt.show()