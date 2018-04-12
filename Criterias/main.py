from scipy.stats import norm
from scipy.stats import kstwobign
import numpy as np
import Smirnov as sm
import LehmanRosenblatt as lr
import AndersonDarling as ad
from Helper import Helper
import xlwt
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt

#если не задавать roundingCount значит без округления
def EmpiricFunctionOfCriteria(n, N, criteria, roundingCount = -1):
    stats = []

    for i in range(N):
        x1 = norm.rvs(loc=0, scale=1, size=n)
        x2 = norm.rvs(loc=0, scale=1, size=n)
        x1_rounded = x1 if roundingCount == -1 else Helper.RoundingArray(x1, roundingCount)
        x2_rounded = x2 if roundingCount == -1 else Helper.RoundingArray(x2, roundingCount)
        stats.append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)

    return ECDF(stats)

def ShowAllEmpiricFunctionsOfSmirnovCriteria(n, m, N, criteria):
    roundingDigitsCounts = [0, 1, 2]  # количество знаков округления значений выборок
    descriptions = ('K-S', '-', '0', '1', '2')
    ecdfs = []
    stats = [[]]
    for digit in roundingDigitsCounts:
        stats.append([])

    for i in range(N):
        x1 = norm.rvs(loc=0, scale=1, size=n)
        x2 = norm.rvs(loc=0, scale=1, size=m)
        stats[0].append(criteria.Result2Samples(x1, x2).statistic)
        for index, digit in enumerate(roundingDigitsCounts):
            x1_rounded = Helper.RoundingArray(x1, digit)
            x2_rounded = Helper.RoundingArray(x2, digit)
            stats[index+1].append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)

    lines = []
    en = np.sqrt(n * m / float(n + m))
    allstats = sum(stats, [])
    x = np.linspace(min(allstats), max(allstats), N)
    lines.append(plt.plot(x, kstwobign.cdf((en + 0.12 + 0.11 / en)*x)))
    for stat in stats:
        ecdf = ECDF(stat)
        ecdfs.append(ecdf)
        lines.append(plt.plot(ecdf.x, ecdf.y))
    plt.legend(labels = descriptions)
    plt.show()

    return ecdfs

def ShowAllEmpiricFunctionsOfLehmRosCriteria(n, m, N, criteria):
    roundingDigitsCounts = [0, 1, 2]  # количество знаков округления значений выборок
    descriptions = ('L-R', '-', '0', '1', '2')
    ecdfs = []
    stats = [[]]
    for digit in roundingDigitsCounts:
        stats.append([])

    for i in range(N):
        x1 = norm.rvs(loc=0, scale=1, size=n)
        x2 = norm.rvs(loc=0, scale=1, size=m)
        stats[0].append(criteria.Result2Samples(x1, x2).statistic)
        for index, digit in enumerate(roundingDigitsCounts):
            x1_rounded = Helper.RoundingArray(x1, digit)
            x2_rounded = Helper.RoundingArray(x2, digit)
            stats[index+1].append(criteria.Result2Samples(x1_rounded, x2_rounded).statistic)

    lines = []
    allstats = sum(stats, [])
    x = np.linspace(min(allstats), max(allstats), N)
    #lines.append(plt.plot(x, lr.LehmanRosenblattCriteria.a1(x)))
    lines.append(plt.plot(x, lr.LehmanRosenblattCriteria.GetStatisticDistribution(x)))
    for stat in stats:
        ecdf = ECDF(stat)
        ecdfs.append(ecdf)
        lines.append(plt.plot(ecdf.x, ecdf.y))
    plt.legend(labels=descriptions)
    plt.show()

    return ecdfs

def GetStatisticsValueByRoundedSeries():
    s_c = sm.SmirnovCriteria()
    lr_c = lr.LehmanRosenblattCriteria()
    ad_c = ad.AndersonDarlingCriteria()

    style0 = xlwt.easyxf('font: name Times New Roman, color-index black, bold on')
    wb = xlwt.Workbook()
    ws = wb.add_sheet('1')
    strNum = 0

    sizes = [10, 20, 30, 50, 100, 200, 500, 1000]#размеры выборок
    roundingDigitsCounts = [0, 1, 2]#количество знаков округления значений выборок
    criterias = {'Смирнов': s_c, 'Леман-Розенблатт': lr_c, 'Андерсон-Дарлинг': ad_c}

    for size in sizes:
        ws.write(strNum, 0, f"n={size}", style0)
        for index, crName in enumerate(criterias):
            ws.write(strNum+index+1, 0, crName)
        x1 = norm.rvs(loc=0, scale=1, size=size)
        x2 = norm.rvs(loc=0, scale=1, size=size)
        for index, digit in enumerate(roundingDigitsCounts):
            x1_rounded = Helper.RoundingArray(x1, digit)
            x2_rounded = Helper.RoundingArray(x2, digit)
            x1_rounded.sort()
            x2_rounded.sort()

            ws.write(strNum, index+1, digit, style0)
            for ind, criteria in enumerate(criterias.values()):
                ws.write(strNum+ind+1, index + 1, criteria.Result2Samples(x1_rounded, x2_rounded).statistic)
        strNum += len(criterias) + 2
    wb.save('Statistics Criterias Rounded.xls')

def SmirnovEmpiricPlot():
    s_c = sm.SmirnovCriteria()
    ecdf = EmpiricFunctionOfCriteria(100, 100, s_c)
    plt.plot(ecdf.x, ecdf.y)
    plt.show()

def SmirnovEmpiricPlots():
    s_c = sm.SmirnovCriteria()
    ecdf = ShowAllEmpiricFunctionsOfSmirnovCriteria(100, 100, 100, s_c)

def LehmRosEmpiricPlots():
    lr_c = lr.LehmanRosenblattCriteria()
    ecdf = ShowAllEmpiricFunctionsOfLehmRosCriteria(100, 100, 100, lr_c)

#GetStatisticsValueByRoundedSeries()
#SmirnovEmpiricPlot()
#SmirnovEmpiricPlots()
LehmRosEmpiricPlots()