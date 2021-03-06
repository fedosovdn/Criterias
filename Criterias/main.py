from scipy.stats import norm
import Smirnov as sm
import LehmanRosenblatt as lr
import AndersonDarling as ad
from Helper import Helper
import xlwt
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from ShowAllEmpiricFunctionsHelper import ShowAllEmpiricFunctionsHelper as shower
from PowerCalculateHelper import PowerCalculateHelper as powCalc
import numpy as np
from scipy.stats import kstwobign

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
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(200, 200, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(500, 500, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(1000, 1000, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(2000, 2000, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(5000, 5000, 16600)

    # x = np.linspace(0.1, 10, 16600)
    # plt.plot(x, kstwobign.cdf(x))
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(200, 200, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(500, 500, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(999, 1001, 16600)
    shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(499, 1501, 16600)
    # shower.ShowAllEmpiricFunctionsOfSmirnovCriteria(5000, 5000, 16600)

    # plt.ylabel("G(Sc|H0)")
    # plt.xlabel("Sc")
    # plt.legend(labels=('K-S', 'n=m=200', 'n=m=500', 'n=m=1000', 'n=m=2000', 'n=m=5000'))
    # plt.show()

def LehmRosEmpiricPlots():
    shower.ShowAllEmpiricFunctionsOfLehmRosCriteria(10000, 10000, 16600, 0, 1)
    # shower.ShowAllEmpiricFunctionsOfLehmRosCriteria(1000, 1000, 16600)
    # shower.ShowAllEmpiricFunctionsOfLehmRosCriteria(2000, 2000, 16600)
    # shower.ShowAllEmpiricFunctionsOfLehmRosCriteria(5000, 5000, 16600)

def ADEmpiricPlots():
    # shower.ShowAllEmpiricFunctionsOfADCriteria(30, 30, 16600)
    # shower.ShowAllEmpiricFunctionsOfADCriteria(30, 40, 16600)
    shower.ShowAllEmpiricFunctionsOfADCriteria(1000, 1000, 16600)
    # shower.ShowAllEmpiricFunctionsOfADCriteria(1000, 1000, 16600)
    # shower.ShowAllEmpiricFunctionsOfADCriteria(2000, 2000, 16600)
    # shower.ShowAllEmpiricFunctionsOfADCriteria(5000, 5000, 16600)

def CalcPowers():
    # cr = ad.AndersonDarlingCriteria()
    cr = lr.LehmanRosenblattCriteria()
    # cr = sm.SmirnovCriteria()
    # powCalc.CalcualteStats(500, 500, 16600, cr, '-')#без округления
    # powCalc.CalcualteStats(2000, 2000, 16600, cr, 1)
    powCalc.CalculateStats(200, 200, 16600, cr, 0)

def ShowDistrPlots():
    powCalc.ShowDistrPlots(500, 0)

def ShowTheorDistrs():
    powCalc.ShowTheorDistrs()

#GetStatisticsValueByRoundedSeries()
#SmirnovEmpiricPlot()
SmirnovEmpiricPlots()
# LehmRosEmpiricPlots()
# ADEmpiricPlots()
# CalcPowers()
# ShowDistrPlots()
# ShowTheorDistrs()