import math


def normalize(point):
    norm = math.sqrt(sum([point[i]**2 for i in range(3)]))
    return tuple([point[i]/norm for i in range(3)])


def locate(point, normalVector, supportVector):
    s = sum([(point[i]-supportVector[i])*normalVector[i] for i in range(3)])
    return s<=0

