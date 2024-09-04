import math
import shapely
from shapely.geometry import Polygon
from matplotlib.patches import Ellipse

def createPointsListForEllipse(x, y, a, b, angle):
    y = -y
    angle = -angle
    theta = [x * 2 * math.pi / 1000 for x in range(1000)]
    coords = []
    for i in theta:
        coords.append([a * math.cos(i) * math.cos(angle) - b * math.sin(i) * math.sin(angle) + x, a * math.cos(i) * math.sin(angle) + b * math.sin(i) * math.cos(angle) + y]) 
    return coords

PEL = Polygon(createPointsListForEllipse(245.9637,354.8639,1139.4186,1374.1214,1.0056))
MXL = Polygon(createPointsListForEllipse(188.7347,343.7469,1298.0515,1327.9438,-1.3532))

total_area = shapely.union_all([PEL, MXL]).area

print(PEL.area)
print(MXL.area)
print(total_area)
print(shapely.intersection_all([PEL, MXL]).area / total_area * 100, "% of variants that are common in one sample are common in all samples.")
print(shapely.symmetric_difference_all([PEL, MXL]).area / total_area * 100, "% of variants are only common in a single sample.")
