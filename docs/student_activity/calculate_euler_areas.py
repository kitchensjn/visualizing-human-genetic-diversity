import math
from shapely.geometry import Polygon
import pandas as pd


def createPointsListForEllipse(x, y, a, b, angle):
    y = -y
    angle = -angle
    theta = [x * 2 * math.pi / 100 for x in range(100)]
    coords = []
    for i in theta:
        coords.append([a * math.cos(i) * math.cos(angle) - b * math.sin(i) * math.sin(angle) + x, a * math.cos(i) * math.sin(angle) + b * math.sin(i) * math.cos(angle) + y]) 
    return coords


def calcIntersectPercentage(data):
    u = Polygon(createPointsListForEllipse(data[0].h, data[0].k, data[0].a, data[0].b, data[0].phi))
    i = Polygon(createPointsListForEllipse(data[0].h, data[0].k, data[0].a, data[0].b, data[0].phi))
    for e in data[1:]:
        ellipse = Polygon(createPointsListForEllipse(e.h, e.k, e.a, e.b, e.phi))
        u = u.union(ellipse)
        i = i.intersection(ellipse)
    return i.area / u.area

def calcUniquePercentage(data):
    all = []
    u = Polygon(createPointsListForEllipse(data[0].h, data[0].k, data[0].a, data[0].b, data[0].phi))
    for i,p in enumerate(data):
        poly_p = Polygon(createPointsListForEllipse(p.h, p.k, p.a, p.b, p.phi))
        if i != 0:
            u = u.union(poly_p)
        for j,q in enumerate(data):
            if i != j:
                poly_q = Polygon(createPointsListForEllipse(q.h, q.k, q.a, q.b, q.phi))
                poly_p = poly_p.difference(poly_q)
        all.append(poly_p.area)
    return sum(all) / u.area



def calculate_areas_of_euler_diagram(row):
    print(row)





guesses = pd.read_csv("guesses.csv")

guesses["ACB.a"] = 1597.628
guesses["ACB.b"] = 1597.628

guesses["CEU.a"] = 1350.0972
guesses["CEU.b"] = 1350.0972

guesses["MXL.a"] = 1342.6299
guesses["MXL.b"] = 1342.6299

guesses.apply(calculate_areas_of_euler_diagram)