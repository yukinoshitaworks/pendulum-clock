"""Fit test coupon: print this before the clock parts to choose ROD_PRESS and BRG_R for your printer and filament.

Row 1 (front): through holes for the phi3 rod, diameters 3.35 / 3.40 / 3.45 / 3.50 / 3.55 mm.
Row 2 (back):  623ZZ pockets 4 mm deep (with the phi6.5 clearance behind, as in the plates), 10.30 .. 10.50 mm.
The diameter is engraved next to each hole. Print it flat, engraved face up, with the same slicer settings as the parts.

    python3 fit_test.py        -> out/extras/fit_test_coupon.step and .stl
"""
import os
import cadquery as cq
from shapely.geometry import Polygon
from shapely import affinity
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

V = cq.Vector
OUT = os.environ.get("CLOCK_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"))
os.makedirs(OUT + "/extras", exist_ok=True)

PRESS = [3.35, 3.40, 3.45, 3.50, 3.55]
POCKET = [10.30, 10.35, 10.40, 10.45, 10.50]
W, H, T = 92.0, 42.0, 5.0
XS = [10.0 + 18.0 * i for i in range(5)]
Y_ROD, Y_BRG = 12.0, 28.0


def cyl(x, y, r, z0, z1):
    return cq.Solid.makeCylinder(r, z1 - z0, V(x, y, z0), V(0, 0, 1))


def text_poly(s, size, cx, cy):
    tp = TextPath((0, 0), s, size=size, prop=FontProperties(family="DejaVu Sans", weight="bold"))
    out = Polygon()
    for p in tp.to_polygons():
        if len(p) >= 3:
            out = out.symmetric_difference(Polygon(p).buffer(0))
    x0, y0, x1, y1 = out.bounds
    return affinity.translate(out, cx - (x0 + x1) / 2, cy - (y0 + y1) / 2)


def engrave(solid, poly, depth):
    polys = list(poly.geoms) if poly.geom_type == "MultiPolygon" else [poly]
    for pg in polys:
        ow = cq.Wire.makePolygon([V(x, y, T - depth) for x, y in list(pg.exterior.coords)[:-1]], close=True)
        iw = [cq.Wire.makePolygon([V(x, y, T - depth) for x, y in list(r.coords)[:-1]], close=True) for r in pg.interiors]
        solid = solid.cut(cq.Solid.extrudeLinear(cq.Face.makeFromWires(ow, iw), V(0, 0, depth + 1)))
    return solid


plate = cq.Solid.makeBox(W, H, T, V(0, 0, 0))
for x, d in zip(XS, PRESS):
    plate = plate.cut(cyl(x, Y_ROD, d / 2, -1, T + 1))
    plate = engrave(plate, text_poly("%.2f" % d, 3.2, x, Y_ROD - 6.5), 0.6)
for x, d in zip(XS, POCKET):
    plate = plate.cut(cyl(x, Y_BRG, d / 2, T - 4.0, T + 1)).cut(cyl(x, Y_BRG, 3.25, -1, T))
    plate = engrave(plate, text_poly("%.2f" % d, 3.2, x, Y_BRG + 9.0), 0.6)

cq.exporters.export(cq.Workplane(obj=plate), OUT + "/extras/fit_test_coupon.step")
cq.exporters.export(cq.Workplane(obj=plate), OUT + "/extras/fit_test_coupon.stl", tolerance=0.01, angularTolerance=0.1)
print("fit_test_coupon: %.0f x %.0f x %.0f mm, %.1f cm3 -> %s/extras/" % (W, H, T, plate.Volume() / 1000, OUT))
