from decimal import *


class Helper:
    @staticmethod
    def GetLastIndexOf(list, value):
        return len(list) - list[::-1].index(value) - 1

    '''Получить ранг элемента вариационного ряда'''
    @staticmethod
    def GetRang(list, value):
        first = list.index(value)
        last = first
        for elem in list[first+1:]:
            if elem == value:
                last += 1
                continue
            break
        return (first + last) / 2

    @staticmethod
    def RoundingArray(array, digitCount):
        # getcontext().rounding = ROUND_HALF_UP
        result = []
        for elem in array:
            result.append(round(Decimal(elem), digitCount))
        return result
