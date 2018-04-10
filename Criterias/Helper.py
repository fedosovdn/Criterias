from decimal import *


class Helper:
    @staticmethod
    def GetLastIndexOf(list, value):
        return len(list) - list[::-1].index(value) - 1

    @staticmethod
    def RoundingArray(array, digitCount):
        # getcontext().rounding = ROUND_HALF_UP
        result = []
        for elem in array:
            #result.append(format(elem, f".{digitCount}f"))
            result.append(round(Decimal(elem), digitCount))
        return result
