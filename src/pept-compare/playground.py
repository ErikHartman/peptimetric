from methods import *
from lists import *


def get_thresholds(lst):
    thresholds = []
    lst.sort()
    for i in range(len(lst)):
        if i % (int(len(lst) / 5)) == 0 and lst[i] != 0 and len(thresholds) <= 5:
            thresholds.append(lst[i])
    return thresholds


lst = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 77, 7, 7, 7
       , 10, 10, 7,8,9,7,5,4,45,6,7,8]
print(len(lst), get_thresholds(lst))
