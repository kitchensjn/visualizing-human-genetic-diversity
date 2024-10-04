import math
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely import intersection_all
import matplotlib.pyplot as plt
from itertools import combinations


def createPointsListForEllipse(x, y, a, b, angle):
    y = -y
    angle = -angle
    theta = [x * 2 * math.pi / 100 for x in range(100)]
    coords = []
    for i in theta:
        coords.append([a * math.cos(i) * math.cos(angle) - b * math.sin(i) * math.sin(angle) + x, a * math.cos(i) * math.sin(angle) + b * math.sin(i) * math.cos(angle) + y]) 
    return coords


def calcIntersectPercentage(data):
    filtered = []
    for ellipse in data:
        if ellipse["active"] == "true":
            filtered.append(ellipse)
    data = filtered
    if len(data) > 1:
        u = Polygon(createPointsListForEllipse(data[0]["h"], data[0]["k"], data[0]["a"], data[0]["b"], data[0]["phi"]))
        i = Polygon(createPointsListForEllipse(data[0]["h"], data[0]["k"], data[0]["a"], data[0]["b"], data[0]["phi"]))
        for e in data[1:]:
            ellipse = Polygon(createPointsListForEllipse(e["h"], e["k"], e["a"], e["b"], e["phi"]))
            u = u.union(ellipse)
            i = i.intersection(ellipse)
        return i.area / u.area
    elif len(data) > 0:
        return 1
    else:
        return 0

def calcUniquePercentage(data):
    filtered = []
    for ellipse in data:
        if ellipse["active"] == "true":
            filtered.append(ellipse)
    data = filtered
    if len(data) > 1:
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
    elif len(data) > 0:
        return 1
    else:
        return 0


data = [
    {"abbreviation":"PEL","h":1,"k":0,"a":1,"b":1,"phi":0,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"white","stroke_dasharray":"none","active":"true"},
    {"abbreviation":"MXL","h":0,"k":1,"a":1,"b":1,"phi":0,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"white","stroke_dasharray":"none","active":"true"},
    {"abbreviation":"CEU","h":0,"k":0.5,"a":1,"b":1,"phi":0,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"white","stroke_dasharray":"none","active":"false"}
]

data = [
    {"abbreviation":"","h":0,"k":-1750,"a":1500,"b":1500,"phi":0,"color":"#E69F00","stroke_dasharray":"none"},
    {"abbreviation":"","h":0,"k":1750,"a":1500,"b":1500,"phi":0,"color":"#56B4E9","stroke_dasharray":"none"}
]


"""
def calculate_intersection_areas(data):
    Calculates the area fraction of all categories (unique to one, unique to two, ...)

    Parameters
    ----------
    data : list
        List of dictionaries which contains the ellipses data

    Returns
    -------
    areas : list
        List of area fractions associated with each category
    

    ellipses = []
    for ellipse in data:
        poly = Polygon(createPointsListForEllipse(ellipse["h"], ellipse["k"], ellipse["a"], ellipse["b"], ellipse["phi"]))
        ellipses.append(poly)
    areas = []
    for a in range(len(ellipses)):
        a_sum = 0
        for combo in combinations(ellipses, a+1):
            poly_combo = intersection_all(combo)
            not_combo = [e for e in ellipses if e not in combo]
            poly_not = unary_union(not_combo)
            diff = poly_combo.difference(poly_not)
            a_sum += diff.area
        areas.append(a_sum)
    total = sum(areas)
    areas = [a/total for a in areas]
    return areas

def new_calculate_intersection_areas(data):
    ellipses = []
    for ellipse in data:
        poly = Polygon(createPointsListForEllipse(ellipse["h"], ellipse["k"], ellipse["a"], ellipse["b"], ellipse["phi"]))
        ellipses.append(poly)
    areas = []
    for a in range(len(ellipses)):
        a_sum = 0
        for combo in combinations(ellipses, a+1):
            poly_combo = combo[0]
            for i in range(len(combo)-1):
                poly_combo = poly_combo.intersection(combo[i+1])
            not_combo = [e for e in ellipses if e not in combo]
            poly_not = unary_union(not_combo)
            diff = poly_combo.difference(poly_not)
            a_sum += diff.area
        areas.append(a_sum)
    total = sum(areas)
    areas = [a/total for a in areas]
    return areas
"""


def createPointsListForEllipse(x, y, a, b, angle):
    y = -y
    angle = -angle
    theta = [x * 2 * math.pi / 100 for x in range(100)]
    coords = []
    for i in theta:
        coords.append([a * math.cos(i) * math.cos(angle) - b * math.sin(i) * math.sin(angle) + x, a * math.cos(i) * math.sin(angle) + b * math.sin(i) * math.cos(angle) + y]) 
    return coords

def calculate_intersection_areas(data):
    ellipses = []
    for ellipse in data:
        poly = Polygon(createPointsListForEllipse(ellipse["h"], ellipse["k"], ellipse["a"], ellipse["b"], ellipse["phi"]))
        ellipses.append(poly)
    areas = []
    for a in range(len(ellipses)):
        a_sum = 0
        for combo in combinations(ellipses, a+1):
            print(combo)
            poly_combo = combo[0]
            for i in range(len(combo)-1):
                poly_combo = poly_combo.intersection(combo[i+1])
            not_combo = [e for e in ellipses if e not in combo]
            poly_not = unary_union(not_combo)
            diff = poly_combo.difference(poly_not)
            print(diff.area)
            a_sum += diff.area
        areas.append(a_sum)
    if len(areas) > 0:
        total = sum(areas)
        areas = [round((a/total)*100) for a in areas]
    else:
        areas = ["~"]
    return areas

print(calculate_intersection_areas(data=data))

for c in combinations([1,1,2], 2):
    print(c)
#print(new_calculate_intersection_areas(data=data))