from enum import Enum

class CriteriaSide(Enum):
    left = 1
    right = 2
    both = 3

class PowerCalculateHelper:
    @staticmethod
    def CalculatePower(statsH0, statsH1, criteriaSide, alphas):
        powers = []

        return powers