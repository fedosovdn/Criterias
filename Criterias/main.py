from scipy.stats import norm
import Smirnov as sm
import LehmanRosenblatt as lr
import AndersonDarling as ad
from Helper import Helper
import xlwt


def GetStatisticsValueByRoundedSeries():
    s_c = sm.SmirnovCriteria()
    lr_c = lr.LehmanRosenblattCriteria()
    ad_c = ad.AndersonDarlingCriteria()

    style0 = xlwt.easyxf('font: name Times New Roman, color-index black, bold on')
    wb = xlwt.Workbook()
    ws = wb.add_sheet('1')
    strNum = 0

    sizes = [50, 100, 200, 500]#размеры выборок
    roundingDigitsCounts = [0, 1, 2]#количество знаков округления значений выборок
    criterias = ["Смирнов","Леман-Розенблатт", "Андерсон-Дарлинг"]

    for size in sizes:
        ws.write(strNum, 0, f"n={size}", style0)
        for index, crName in enumerate(criterias):
            ws.write(strNum+index+1, 0, crName)
        x1 = norm.rvs(loc=0, scale=1, size=size)
        x2 = norm.rvs(loc=0, scale=1, size=size)
        for index, digit in enumerate(roundingDigitsCounts):
            x1_rounded = Helper.RoundingArray(x1, digit)
            x2_rounded = Helper.RoundingArray(x2, digit)
            x1_rounded.sort()#важно! должны быть сортированы
            x2_rounded.sort()

            ws.write(strNum, index+1, digit, style0)
            ws.write(strNum+1, index + 1, s_c.SciPyResult2Samples(x1_rounded, x2_rounded).statistic)
            ws.write(strNum+2, index + 1, lr_c.Result(x1_rounded, x2_rounded))
            ws.write(strNum+3, index + 1, ad_c.SciPiResultKSamples([x1_rounded, x2_rounded]).statistic)
        strNum += len(criterias) + 2
    wb.save('Statistics Criterias Rounded.xls')


GetStatisticsValueByRoundedSeries()