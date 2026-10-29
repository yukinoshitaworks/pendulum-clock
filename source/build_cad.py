import os, sys, json, math, struct, base64, itertools
from math import pi, sin, cos, tan, radians, degrees, sqrt, atan2
import numpy as np
import cadquery as cq
from shapely.geometry import Polygon, Point, LineString, box as sbox
from shapely.ops import unary_union
from shapely import affinity
import design as D

V = cq.Vector
OUT = os.environ.get("CLOCK_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"))
os.makedirs(OUT + "/parts_step", exist_ok=True)
os.makedirs(OUT + "/parts_stl", exist_ok=True)

# ------------------------------------------------------------------ helpers
def cyl(x, y, r, z0, z1):
    return cq.Solid.makeCylinder(r, z1 - z0, V(x, y, z0), V(0, 0, 1))

def cyl_dir(p, d, r, h):
    return cq.Solid.makeCylinder(r, h, V(*p), V(*d))

def tube(x, y, ro, ri, z0, z1):
    return cyl(x, y, ro, z0, z1).cut(cyl(x, y, ri, z0 - 1, z1 + 1))

def cone(x, y, r0, r1, z0, z1):
    return cq.Solid.makeCone(r0, r1, z1 - z0, V(x, y, z0), V(0, 0, 1))

def prism_poly(poly, z0, z1, at=(0, 0)):
    poly = poly.simplify(0.02)
    polys = list(poly.geoms) if poly.geom_type == "MultiPolygon" else [poly]
    out = None
    for pg in polys:
        ow = cq.Wire.makePolygon([V(x + at[0], y + at[1], z0) for x, y in list(pg.exterior.coords)[:-1]], close=True)
        iw = [cq.Wire.makePolygon([V(x + at[0], y + at[1], z0) for x, y in list(r.coords)[:-1]], close=True) for r in pg.interiors]
        s = cq.Solid.extrudeLinear(cq.Face.makeFromWires(ow, iw), V(0, 0, z1 - z0))
        out = s if out is None else out.fuse(s)
    return out

def prism_segs(segs, z0, z1, at=(0, 0)):
    edges = []
    v = lambda p: V(p[0] + at[0], p[1] + at[1], z0)
    for kind, p in segs:
        if kind == "line":
            edges.append(cq.Edge.makeLine(v(p[0]), v(p[1])))
        elif kind == "arc":
            edges.append(cq.Edge.makeThreePointArc(v(p[0]), v(p[1]), v(p[2])))
        else:
            edges.append(cq.Edge.makeSpline([v(q) for q in p]))
    w = cq.Wire.assembleEdges(edges)
    return cq.Solid.extrudeLinear(cq.Face.makeFromWires(w, []), V(0, 0, z1 - z0))

def gear_solid(g, z0=None, z1=None):
    return prism_segs(g.segments(), g.z0 if z0 is None else z0, g.z1 if z1 is None else z1, at=g.center)

def hexagon(af, rot=0.0):
    r = af / 2 / cos(pi / 6)
    return Polygon([(r * cos(rot + i * pi / 3), r * sin(rot + i * pi / 3)) for i in range(6)])

def sector(r0, r1, a0, a1, n=24):
    o = [(r1 * cos(a0 + (a1 - a0) * i / n), r1 * sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]
    i_ = [(r0 * cos(a1 - (a1 - a0) * i / n), r0 * sin(a1 - (a1 - a0) * i / n)) for i in range(n + 1)]
    return Polygon(o + i_)

def spoke_pockets(center, r_in, r_out, n, spoke_w, z0, z1, phase=0.0, fillet=1.6):
    """n lightening pockets leaving n straight spokes of width spoke_w."""
    ring = Point(0, 0).buffer(r_out, 24).difference(Point(0, 0).buffer(r_in, 16))
    for k in range(n):
        a = phase + k * 2 * pi / n
        sp = LineString([(0, 0), ((r_out + 5) * cos(a), (r_out + 5) * sin(a))]).buffer(spoke_w / 2, 2, cap_style=2)
        ring = ring.difference(sp)
    ring = ring.buffer(-fillet, 4).buffer(fillet, 4)
    return prism_poly(ring, z0 - 1, z1 + 1, at=center)

def bars(segs, w, extra=()):
    shp = [LineString([a, b]).buffer(w / 2, 8) for a, b in segs] + [Point(c).buffer(r, 10) for c, r in extra]
    u = unary_union(shp)
    return u.buffer(3.0, 5).buffer(-3.0, 5)

def clean(s):
    try:
        return s.clean()
    except Exception:
        return s

# ------------------------------------------------------------------ registry
PARTS = []   # dict(name, ja, solid, group, color, kind, qty_of, print_rot)
def add(name, ja, solid, group="static", color="plate", kind="print", print_rot=None, share=None, note=""):
    PARTS.append(dict(name=name, ja=ja, solid=clean(solid), group=group, color=color, kind=kind,
                      print_rot=print_rot, share=share, note=note))

C, G, E, T, A, M, H, P, D0 = D.C, D.G, D.E, D.T, D.A, D.M, D.H, D.P, D.D0
GAP, PT = D.GAP, D.PLATE_T
BELL = D.BELL
# Hole sizes (diameters). 2026-10-12: set after the first prints came out too tight. 2026-10-13: nut and key hexes opened up. 2026-10-17: rod press 3.3, bearing pocket 10.3, running hole 3.6.
# 2026-10-21: test assembly on a Creality K1 with Hyper PLA-CF: press holes and bearing pockets still too tight.
# K1 prints small holes 0.3-0.7 mm under size (more for smaller holes), so press 3.45 and pocket 10.4.
# 2026-10-23: 3.45 turned out too loose for a press fit; 3.40.
ROD_PRESS = 3.40 / 2   # phi3 rod pressed into a part, plate or the pendulum hanger (was 2.95 / 3.0; 3.2 on 10-12, 3.3 on 10-17, 3.45 on 10-21)
ROD_FREE  = 3.60 / 2   # part that turns on a phi3 rod: great wheel, drum, minute wheel, hammer (was 3.3; 3.5 on 10-12)
BRG_R     = 10.4 / 2   # 623ZZ pocket (was 10.1; 10.2 on 10-12, 10.3 on 10-17)
M3_R      = 3.50 / 2   # M3 clearance: plates, pillars, standoffs, spacers, wall frame, dial, drum flange (was 3.4)
TAP3_R    = 2.80 / 2   # M3 self-tapping pilot (was 2.6)
TAP4_R    = 3.50 / 2   # M4 pilot: bell post, pendulum rod into the hanger (was 3.3 / 3.4)
X         = 0.2 / 2    # every other bore and pin hole: +0.2 on the diameter
PEND_M4 = 0.2 / 2      # 2026-10-23: M4 holes of the pendulum (rod into the hanger, bob and rating-nut bores) another +0.2
HANG_ROD_R = 4.5 / 2   # 2026-10-30: rod hole in the hanger phi 4.5 (3.7 on 10-23, 3.9 on 10-26), as asked
DRUM_TAP_R = 3.10 / 2  # 2026-10-28: M3 tapping screws into the drum were too tight (stripped while tightening): 2.8 -> 3.1
HEX_M3, HEX_M4, HEX_KEY = 5.9, 7.3, 10.5   # across flats: M3 nut (wall frame), M4 nut (rating nut), winding key (was 5.8 / 7.2 / 10.4)
TIP_R = 2.7

# ------------------------------------------------------------------ frame
frame_segs = [((0, G[1]), (0, 116)), (P[0], P[1]), (T, E), (H, E), (C, P[2]), (C, P[3]),
              (P[0], T), (T, P[2]), (P[1], H), (H, P[3])]
bosses = [(C, 9.5), (T, 9.5), (E, 9.5), (A, 9.5), (H, 9.5), (G, 12.0)] + [(p, 7.5) for p in P]

back2d = unary_union([bars(frame_segs + [(H, BELL), (P[3], BELL)], 13, bosses + [(BELL, 8.5)])])
back = prism_poly(back2d, -PT, 0)
for c in (C, T, E, A):
    back = back.cut(cyl(c[0], c[1], BRG_R, -4, 0.1)).cut(cyl(c[0], c[1], 3.25, -PT - 1, 0))
back = back.cut(cyl(G[0], G[1], ROD_PRESS, -PT - 1, 1))
for p in P:
    back = back.cut(cyl(p[0], p[1], M3_R, -PT - 1, 1))
BELL_RIM_Z, BELL_H, BELL_TH = 18.0, 21.0, 1.5
BELL_RS = (D.BELL_R**2 + BELL_H**2) / (2 * BELL_H)
POST_TOP = BELL_RIM_Z + BELL_H - BELL_TH - 0.7
back = back.fuse(cyl(BELL[0], BELL[1], 6.0, -0.5, POST_TOP)).cut(cyl(BELL[0], BELL[1], TAP4_R, POST_TOP - 14, POST_TOP + 1))
add("back_plate", "後ろ地板（ベル支柱つき）", back, color="plate", print_rot=None)

front2d = bars(frame_segs, 13, bosses + [(M, 6.5)])
front = prism_poly(front2d, GAP, GAP + PT)
for c in (C, T, E, A):
    front = front.cut(cyl(c[0], c[1], BRG_R, GAP - 0.1, GAP + 4)).cut(cyl(c[0], c[1], 3.25, GAP, GAP + PT + 1))
front = front.cut(cyl(G[0], G[1], 6.2 + X, GAP - 1, GAP + PT + 1))
for c in (H, M):
    front = front.cut(cyl(c[0], c[1], ROD_PRESS, GAP - 1, GAP + PT + 1))
front = front.cut(cyl(D0[0], D0[1], M3_R, GAP - 1, GAP + PT + 1)).cut(cone(D0[0], D0[1], 3.3, M3_R, GAP - 0.01, GAP + 1.6))
for p in P:
    front = front.cut(cyl(p[0], p[1], M3_R, GAP - 1, GAP + PT + 1))
add("front_plate", "前地板", front, color="plate", print_rot=("x", 180))

REAR_L = 30.0
Z_WALL0, Z_WALL1 = -PT - REAR_L - 6.0, -PT - REAR_L
TIE_Z0 = Z_WALL0 + 0.6        # rear end of the M3 threaded rods (nut pocket bottom face is at Z_WALL0 + 3)
ZO_ = GAP + PT
for i, p in enumerate(P):
    add(f"pillar_{i+1}", "地板間ピラー (36 mm)", tube(p[0], p[1], 5, M3_R, 0, GAP), color="pillar", share="pillar_36_x4")
    add(f"rear_standoff_{i+1}", "背面スタンドオフ (30 mm)", tube(p[0], p[1], 5, M3_R, -PT - REAR_L, -PT), color="pillar", share="rear_standoff_30_x4")
    # M3 threaded rod: stock 100 mm for all four (rear end in the wall-frame nut pocket, front nut fully engaged)
    add(f"tie_rod_{i+1}", "M3 全ねじ 100 mm", cyl(p[0], p[1], 1.5, TIE_Z0, TIE_Z0 + 100.0), color="steel", kind="hardware")
    if i < 2:   # top two: a printed sleeve takes up the extra length so the front nut sits at the rod end
        add(f"front_spacer_{i+1}", "上側の前スペーサ (%.1f mm)" % (TIE_Z0 + 100.0 - 2.4 - ZO_), tube(p[0], p[1], 5, M3_R, ZO_, TIE_Z0 + 100.0 - 2.4), color="pillar", share="front_spacer_16_x2")

WT = [(-22.0, 131.0), (22.0, 131.0)]
wall2d = bars([(P[0], P[1]), (P[2], P[3]), (P[0], P[2]), (P[1], P[3]), (P[0], WT[0]), (P[1], WT[1]), (WT[0], WT[1])], 14,
              [(p, 7.5) for p in P] + [(w, 8.0) for w in WT])
wall = prism_poly(wall2d, Z_WALL0, Z_WALL1)
for p in P:
    wall = wall.cut(cyl(p[0], p[1], M3_R, Z_WALL0 - 1, Z_WALL1 + 1)).cut(prism_poly(hexagon(HEX_M3), Z_WALL0 - 1, Z_WALL0 + 3.0, at=p))
for w in WT:
    wall = wall.cut(cyl(w[0], w[1], 2.25 + X, Z_WALL0 - 1, Z_WALL1 + 1))
add("wall_frame", "壁掛けフレーム", wall, color="plate", print_rot=("x", 180))

# ------------------------------------------------------------------ bearings & arbors (hardware)
def bearing(c, z0):
    return tube(c[0], c[1], 5.0, 1.5, z0, z0 + 4)
for nm, c in (("C", C), ("T", T), ("E", E), ("A", A)):
    add(f"bearing_{nm}_rear", "ベアリング 623ZZ", bearing(c, -4), color="steel", kind="hardware")
    add(f"bearing_{nm}_front", "ベアリング 623ZZ", bearing(c, GAP), color="steel", kind="hardware")
# phi 3 rods at stock lengths (5 mm steps): 70 / 45 / 45 / 65 / 55 / 20 / 20
RODS = {"C": (C, -4.0, 66.0, "centre"), "T": (T, -5.0, 40.0, "third"), "E": (E, -5.0, 40.0, "escape"),
        "A": (A, -25.0, 40.0, "anchor"), "G": (G, -5.0, 50.0, "static"), "M": (M, GAP, 56.0, "static"), "H": (H, 21.0, GAP + PT, "static")}
for nm, (c, z0, z1, grp) in RODS.items():
    add(f"rod_{nm}", f"φ3 丸棒 {z1-z0:.1f} mm", cyl(c[0], c[1], 1.5, z0, z1), group=grp, color="steel", kind="hardware")

# ------------------------------------------------------------------ great wheel, click, drum
gw = gear_solid(D.great).fuse(cyl(G[0], G[1], 6.0, 0.4, 9.0))
gw = gw.cut(spoke_pockets(G, 32.0, D.great.rf - 4.0, 6, 6.0, 9, 14, phase=radians(30)))
PEG_R, PEG_A = 26.0, [radians(60 + 120 * k) for k in range(3)]
for a in PEG_A:
    gw = gw.cut(cyl(G[0] + PEG_R * cos(a), G[1] + PEG_R * sin(a), 2.6 + X, 8, 15))
gw = gw.cut(cyl(G[0], G[1], ROD_FREE, -1, 20))
add("great_wheel", "一番車 60枚 (m1.5)", gw, group="great", color="gear1", print_rot=("x", 180))

# 2026-10-23: the ratchet slipped backwards under the weight. The tooth faces were radial, so the load had no
# component holding the pawls in and only the thin spring arms kept them engaged. Now the faces are undercut by 12 deg
# (load pulls the pawl into the tooth), the teeth are 3.0 mm deep (was 2.5), the pawls follow the undercut face and the
# spring arms are 1.8 mm wide (was 1.6) and start from the back of the pawl, clear of the deeper teeth.
RAT_N, RAT_TIP, RAT_ROOT = 12, 15.5, 12.5
RAT_PH = radians(10)
RAT_HOOK_DEG = 12.0
def _hook_angle():
    # angular offset of the tooth tip so that the face (root -> tip) leans RAT_HOOK_DEG from radial, over the pawl
    lo, hi = 0.0, 0.3
    for _ in range(60):
        b = (lo + hi) / 2
        tx, ty = RAT_TIP * cos(b) - RAT_ROOT, RAT_TIP * sin(b)
        lo, hi = (b, hi) if degrees(atan2(ty, tx)) < RAT_HOOK_DEG else (lo, b)
    return lo
RAT_HOOK = _hook_angle()
def ratchet2d():
    pts = []
    for k in range(RAT_N):
        a = RAT_PH + k * 2 * pi / RAT_N
        pts += [(RAT_TIP * cos(a + RAT_HOOK), RAT_TIP * sin(a + RAT_HOOK)), (RAT_ROOT * cos(a), RAT_ROOT * sin(a))]
    return Polygon(pts)

def click2d():
    parts = [Point(0, 0).buffer(29.5, 24).difference(Point(0, 0).buffer(22.5, 20))]
    pol = lambda r, a: (r * cos(a), r * sin(a))
    rot = lambda p, a: (p[0] * cos(a) - p[1] * sin(a), p[0] * sin(a) + p[1] * cos(a))
    # point on the tooth face line (in the frame of a face at angle 0) at radius r
    fx, fy = RAT_TIP * cos(RAT_HOOK) - RAT_ROOT, RAT_TIP * sin(RAT_HOOK); fl = sqrt(fx * fx + fy * fy); ux, uy = fx / fl, fy / fl
    def on_face(r):
        bq = 2 * RAT_ROOT * ux; cq_ = RAT_ROOT ** 2 - r * r; t = (-bq + sqrt(bq * bq - 4 * cq_)) / 2
        return (RAT_ROOT + t * ux, t * uy)
    gap = radians(0.6)
    for k in range(3):
        f = RAT_PH + k * 2 * pi / 3
        arm = [pol(16.2 + (23.0 - 16.2) * i / 26, f + radians(9 + 99 * i / 26)) for i in range(27)]
        parts.append(LineString(arm).buffer(0.9, 3))
        parts.append(Polygon([rot(on_face(12.9), f + gap), rot(on_face(15.9), f + gap), pol(16.6, f + radians(11)), pol(14.35, f + radians(14))]))
    return unary_union(parts)

Z_CL0, Z_CL1 = 14.3, 19.6
click = prism_poly(click2d(), Z_CL0, Z_CL1, at=G)
for a in PEG_A:
    click = click.fuse(cyl(G[0] + PEG_R * cos(a), G[1] + PEG_R * sin(a), 2.45, 9.4, Z_CL0 + 0.5))
add("click_plate", "コハゼ板（一体ばね爪×3）", click, group="great", color="accent", print_rot=("x", 180))

# 2026-10-26: barrel r 20.0 -> 23.5 (cord radius 24.0), flanges r 23.5 -> 27.5. The cone under the rear flange keeps
# about the same overhang as before (4 mm over 3 mm); the barrel is 7.0 mm wide (was 7.5).
FL_R = D.DRUM_R + 4.0
drum_a = prism_poly(ratchet2d(), 14.6, 20.1, at=G).fuse(cyl(G[0], G[1], FL_R, 20.0, 21.5)) \
    .fuse(cone(G[0], G[1], FL_R, D.DRUM_R, 21.5, 24.5)).fuse(cyl(G[0], G[1], D.DRUM_R, 24.5, 31.5))
drum_b = cyl(G[0], G[1], FL_R, 31.5, 34.0).fuse(cyl(G[0], G[1], 8.0, 34.0, 35.6)).fuse(cyl(G[0], G[1], 6.0, 35.6, 41.6)) \
    .fuse(prism_poly(hexagon(10.0), 41.6, 53.0, at=G))
for k in range(3):
    a = radians(90 + 120 * k)
    x, y = G[0] + 12 * cos(a), G[1] + 12 * sin(a)
    drum_a = drum_a.cut(cyl(x, y, DRUM_TAP_R, 22.0, 32))
    drum_b = drum_b.cut(cyl(x, y, M3_R, 31, 35)).cut(cone(x, y, M3_R, 3.3, 32.4, 34.01))
CORD_HOLE_R = D.DRUM_R + 1.6   # cord tie hole just outside the barrel
drum_b = drum_b.cut(cyl(G[0] + CORD_HOLE_R * cos(radians(200)), G[1] + CORD_HOLE_R * sin(radians(200)), 1.0, 31, 35))
drum_a = drum_a.cut(cyl(G[0], G[1], ROD_FREE, 10, 40))
drum_b = drum_b.cut(cyl(G[0], G[1], ROD_FREE, 30, 60))
add("drum_core", "巻胴（ラチェットつき）", drum_a, group="drum", color="drum", print_rot=("x", 180))
add("drum_front", "巻胴 前フランジ（六角巻き軸）", drum_b, group="drum", color="drum")

# ------------------------------------------------------------------ going train
# pockets start outside the centre pinion so the pinion (21 teeth since 2026-10-23) never overhangs a pocket
cp = gear_solid(D.c_whl).cut(spoke_pockets(C, max(14.5, D.c_pin.ra + 1.5), D.c_whl.rf - 3.5, 6, 5.0, 2, 7))
cp = cp.fuse(gear_solid(D.c_pin)).fuse(cyl(C[0], C[1], 4.0, 15.0, 23.8)).cut(cyl(C[0], C[1], ROD_PRESS, 0, 30))
add("centre_wheel", "二番車 96枚＋カナ24枚", cp, group="centre", color="gear2")
add("spacer_C_rear", "スペーサ 1.7 mm", tube(C[0], C[1], TIP_R, 1.6 + X, 0.3, 2.0), group="centre", color="pillar")

cam = prism_poly(D.cam_polygon(), D.CAM_Z[0], D.CAM_Z[1], at=C).fuse(cyl(C[0], C[1], 4.0, D.CAM_Z[1], 35.7)).cut(cyl(C[0], C[1], ROD_PRESS, 20, 40))
add("strike_cam", "打鐘カム（1時間で1回転の渦巻き）", cam, group="centre", color="accent")

tp = gear_solid(D.t_whl).cut(spoke_pockets(T, 9.5, D.t_whl.rf - 3.5, 5, 4.5, 15.5, 20.5, phase=radians(18)))
tp = tp.fuse(gear_solid(D.t_pin)).fuse(cyl(T[0], T[1], TIP_R, 0.3, 1.6)).cut(cyl(T[0], T[1], ROD_PRESS, -1, 30))
add("third_wheel", "三番車 80枚＋カナ8枚", tp, group="third", color="gear3", print_rot=("x", 180))
add("spacer_T_front", "スペーサ 15.2 mm", tube(T[0], T[1], TIP_R, 1.6 + X, 20.5, 35.7), group="third", color="pillar")

ep = prism_poly(D.esc_polygon(), 9.0, 14.0, at=E).cut(spoke_pockets(E, 8.5, 20.5, 5, 4.0, 9, 14, phase=radians(90)))
ep = ep.fuse(gear_solid(D.e_pin)).fuse(cyl(E[0], E[1], TIP_R, 21.9, 35.7)).cut(cyl(E[0], E[1], ROD_PRESS, 0, 40))
add("escape_wheel", "ガンギ車 30枚＋カナ8枚", ep, group="escape", color="escape")
add("spacer_E_rear", "スペーサ 8.7 mm", tube(E[0], E[1], TIP_R, 1.6 + X, 0.3, 9.0), group="escape", color="pillar")

# 2026-10-23: the 27 mm press along the whole anchor bore was too tight. Press only in the pallet body (z 8.5-14.5);
# the hub towards the front bearing gets the running size as a relief.
def anchor_solid(poly):
    return prism_poly(poly, 8.5, 14.5, at=A).fuse(cyl(A[0], A[1], 4.0, 14.4, 35.7)).cut(cyl(A[0], A[1], ROD_PRESS, 0, 14.5)).cut(cyl(A[0], A[1], ROD_FREE, 14.5, 40))
an = anchor_solid(D.anchor_polygon())
# 2026-10-28: spare anchors for adjusting the depth of engagement on the real clock (extras/, not in the parts list).
# "deep_0.2" is laid out for a wheel 0.2 mm further away, so it reaches 0.2 mm (centre distance) deeper into the
# real wheel; use it if the escape wheel runs through (teeth slip past the pallets). "shallow_0.2" if it jams.
ANCHOR_VARIANTS = {"deep_0.2": 0.2, "deep_0.4": 0.4, "shallow_0.2": -0.2}
add("anchor", "アンクル（グラハム式デッドビート）", an, group="anchor", color="anchor")
add("spacer_A_rear", "スペーサ 8.2 mm", tube(A[0], A[1], TIP_R, 1.6 + X, 0.3, 8.5), group="anchor", color="pillar")

# ------------------------------------------------------------------ motion work & hands
ZO = D.Z_FRONT_OUT
can = gear_solid(D.cannon).fuse(cyl(C[0], C[1], 3.0, 46.5, 63.5))
can = can.cut(cq.Solid.makeBox(8, 4, 4, V(-4, 2.4, 59.9))).cut(cyl(C[0], C[1], ROD_PRESS, 40, 70))
add("cannon_pinion", "筒カナ 12枚（摩擦で軸に固定）", can, group="cannon", color="gear4")
mw = gear_solid(D.m_whl).fuse(gear_solid(D.m_pin)).cut(cyl(M[0], M[1], ROD_FREE, 40, 60))
for k in range(4):
    a = radians(45 + 90 * k)
    mw = mw.cut(cyl(M[0] + 11.3 * cos(a), M[1] + 11.3 * sin(a), 3.4, 41, 46.001))
add("minute_wheel", "日の裏車 36枚＋カナ10枚", mw, group="minwheel", color="gear3")
add("collar_M", "抜け止めカラー", tube(M[0], M[1], 3.5, ROD_PRESS, 53.3, 56.0), color="pillar", share="collar_x2")
hw = gear_solid(D.h_whl).fuse(cyl(C[0], C[1], 4.5, 51.9, 59.5)).cut(cyl(C[0], C[1], 3.2 + X, 45, 62))
for k in range(5):
    a = radians(90 + 72 * k)
    hw = hw.cut(cyl(C[0] + 11.8 * cos(a), C[1] + 11.8 * sin(a), 3.6, 47, 52.001))
add("hour_wheel", "筒車 40枚（時針）", hw, group="hour", color="gear2")

def hand2d(length, w0, w1, tail, hub, tip_style):
    blade = Polygon([(-w0 / 2, 0), (w0 / 2, 0), (w1 / 2, length - 3), (0, length), (-w1 / 2, length - 3)])
    parts = [blade, Point(0, 0).buffer(hub, 12), Polygon([(-w0 / 2, 0), (w0 / 2, 0), (w0 * 0.35, -tail), (-w0 * 0.35, -tail)]),
             Point(0, -tail).buffer(w0 * 0.75, 8)]
    if tip_style == "spade":
        parts.append(Point(0, length * 0.72).buffer(5.2, 10))
    return unary_union(parts)

mh2 = hand2d(78, 5.0, 1.8, 16, 6.5, "plain")
dhole = Point(0, 0).buffer(3.08 + X, 12).intersection(sbox(-5, -5, 5, 2.5 + X))
mh = prism_poly(mh2.difference(dhole), 60.5, 62.9, at=C)
add("minute_hand", "長針", mh, group="cannon", color="hand")
hh2 = hand2d(50, 6.5, 2.6, 12, 8.0, "spade").difference(Point(0, 0).buffer(4.5 + X, 16)).difference(Point(0, 50 * 0.72).buffer(2.6, 8))
add("hour_hand", "短針", prism_poly(hh2, 56.5, 58.9, at=C), group="hour", color="hand")

dial = tube(0, 0, 88, 70, 53.0, 56.0)
for p in (P[2], P[3]):
    dial = dial.fuse(cyl(p[0], p[1], 7.5, 53.0, 56.0))
marks = []
for k in range(12):
    a = pi / 2 - k * pi / 6
    q = k % 3 == 0
    r0, r1, w = (73.0, 85.5, 4.2) if q else (77.0, 85.5, 2.2)
    if k == 0:
        for off in (-3.2, 3.2):
            marks.append(affinity.rotate(sbox(r0, -1.5 + off, r1, 1.5 + off), a, origin=(0, 0), use_radians=True))
    else:
        marks.append(affinity.rotate(sbox(r0, -w / 2, r1, w / 2), a, origin=(0, 0), use_radians=True))
dial = dial.fuse(prism_poly(unary_union(marks), 55.9, 57.0))
for p in (P[2], P[3], D0):
    dial = dial.cut(cyl(p[0], p[1], M3_R, 52, 58))
dial = dial.cut(cone(D0[0], D0[1], M3_R, 3.3, 55.4, 57.01))
# 2026-10-21: winding hole. The key boss (r 9.5) on the drum hex at G cut 2.5 mm into the ring's inner edge.
KEY_R = 9.5
dial = dial.cut(cyl(G[0], G[1], KEY_R + 1.5, 52.0, 58.0))
add("dial_ring", "文字盤リング", dial, color="dial")
for i, p in enumerate((P[2], P[3])):
    add(f"dial_standoff_{i+1}", "文字盤スタンドオフ (12 mm)", tube(p[0], p[1], 5, M3_R, ZO, 53.0), color="pillar", share="dial_standoff_12_x2")
add("dial_post_top", "文字盤ポスト 12時側 (下穴φ%.1f)" % (2 * TAP3_R), tube(D0[0], D0[1], 4.0, TAP3_R, ZO, 53.0), color="pillar")

# ------------------------------------------------------------------ strike
hm = prism_poly(D.hammer_polygon(), 24.0, 29.0, at=H).fuse(cyl(H[0], H[1], 5.0, 28.9, 35.6)) \
    .fuse(cyl(D.HEAD_C[0], D.HEAD_C[1], D.HEAD_R, 24.0, 34.0)).cut(cyl(H[0], H[1], ROD_FREE, 20, 40)).cut(cyl(D.HEAD_C[0], D.HEAD_C[1], 3.15 + X, 20, 40))
add("hammer", "ハンマー（尾がカムに乗る）", hm, group="hammer", color="hammer")
add("collar_H", "抜け止めカラー", tube(H[0], H[1], 3.5, ROD_PRESS, 21.2, 23.8), color="pillar", share="collar_x2")
zc = BELL_RIM_Z + BELL_H - BELL_RS
bell = cq.Solid.makeSphere(BELL_RS, V(BELL[0], BELL[1], zc), V(0, 0, 1), -90, 90, 360) \
    .cut(cq.Solid.makeSphere(BELL_RS - BELL_TH, V(BELL[0], BELL[1], zc), V(0, 0, 1), -90, 90, 360)) \
    .cut(cq.Solid.makeBox(200, 200, 200, V(BELL[0] - 100, BELL[1] - 100, BELL_RIM_Z - 200))).cut(cyl(BELL[0], BELL[1], 2.25, 30, 45))
add("bell_dome", "ベル（金属ドーム φ60・市販品）", bell, group="bell", color="bell", kind="reference")

# ------------------------------------------------------------------ pendulum
PZ = D.PEND_Z
hz0, hz1 = -20.5, -7.0
hang = cq.Solid.makeBox(16, 41, hz1 - hz0, V(-8, A[1] - 30, hz0))
hang = hang.cut(cyl(A[0], A[1], ROD_PRESS, hz0 - 1, hz1 + 1)).cut(cq.Solid.makeBox(1.2, 12, 20, V(-0.6, A[1], hz0 - 1)))
hang = hang.cut(cyl_dir((-10, A[1] + 6.5, PZ), (1, 0, 0), 1.65 + X, 20)).cut(cyl_dir((0, A[1] - 31, PZ), (0, 1, 0), HANG_ROD_R, 21))

# 2026-10-17: the hanger and the bob are printed as two halves split at z = PZ, the plane that holds the axes of their
# horizontal holes. Each half prints split face up, so the holes become open half-round grooves (no support inside),
# and the halves are glued. Two pegs on the rear half locate the front half.
PEG_R, PEG_H = 1.5, 2.0
# 2026-10-22: the peg hole was PEG_R + X (phi 3.2) and the printed peg did not go in on the K1. A printed phi-3 peg in a
# printed hole behaves like the phi-3 rod in a running hole, which was fine at 3.6, so use the same size.
PEG_HOLE_R = ROD_FREE
def split_glue(solid, z, pegs, ja, name, **kw):
    big = 500.0
    rear = solid.intersect(cq.Solid.makeBox(big, big, big, V(-big / 2, -big / 2, z - big)))
    front = solid.intersect(cq.Solid.makeBox(big, big, big, V(-big / 2, -big / 2, z)))
    for x, y in pegs:
        rear = rear.fuse(cyl(x, y, PEG_R, z - 0.5, z + PEG_H))
        front = front.cut(cyl(x, y, PEG_HOLE_R, z - 0.1, z + PEG_H + 0.4))
    add(name + "_rear", ja + " 後ろ半分（ピンつき・接着）", rear, **kw)
    add(name + "_front", ja + " 前半分（接着）", front, print_rot=("x", 180), **kw)

split_glue(hang, PZ, [(5.2, A[1] - 22.0), (-4.8, A[1] - 4.0)],   # 2026-10-30: lower peg 4.8 -> 5.2 off the wider rod groove
           "振り子ハンガー", "pendulum_hanger", group="anchor", color="anchor")
add("pendulum_rod", "M4 全ねじ 285 mm（ハンガーに15 mm差し込む）", cyl_dir((0, A[1] - 300, PZ), (0, 1, 0), 2.0, 285), group="anchor", color="steel", kind="hardware")
BOB_C = (0.0, A[1] - D.BOB_Y)
bob = cyl(BOB_C[0], BOB_C[1], 35.0, PZ - 8, PZ + 8).cut(cyl_dir((0, BOB_C[1] - 40, PZ), (0, 1, 0), 2.3 + X + PEND_M4, 80))
for sx in (-1, 1):
    bob = bob.cut(cyl(sx * 18.0, BOB_C[1], 12.3, PZ - 4, PZ + 9))
split_glue(bob, PZ, [(18.0, BOB_C[1] + 22.0), (-18.0, BOB_C[1] - 22.0)], "振り玉 φ70", "bob", group="anchor", color="bob")
ny1 = BOB_C[1] - 35.0 - 0.5
nut = cyl_dir((0, ny1 - 9, PZ), (0, 1, 0), 9.0, 9).cut(cyl_dir((0, ny1 - 10, PZ), (0, 1, 0), 2.2 + X + PEND_M4, 12))
hexp = prism_poly(hexagon(HEX_M4), 0, 3.6).rotate(V(0, 0, 0), V(1, 0, 0), -90).translate(V(0, ny1 - 3.5, PZ))
add("rating_nut", "歩度調整ナット（M4ナット内蔵）", nut.cut(hexp), group="anchor", color="accent", print_rot=("x", 90))

# ------------------------------------------------------------------ weight & key
WX, WZ = D.WEIGHT_X, D.WEIGHT_Z
W_TOP, W_LEN, W_R = -128.0, 200.0, 25.0
wsh = cyl_dir((WX, W_TOP - W_LEN, WZ), (0, 1, 0), W_R, W_LEN).cut(cyl_dir((WX, W_TOP - W_LEN + 3, WZ), (0, 1, 0), W_R - 2.0, W_LEN))
wsh = wsh.cut(cyl(WX, W_TOP - 8, 1.65 + X, WZ - 30, WZ + 30))
add("weight_shell", "おもり容器 φ50×200（鉄球などを充填）", wsh, group="weight", color="weight", print_rot=("x", 90))
add("weight_pin", "φ3 丸棒 55 mm（ひも掛け）", cyl(WX, W_TOP - 8, 1.5, WZ - 27.5, WZ + 27.5), group="weight", color="steel", kind="hardware")
KX, KY = -125.0, -60.0
# 2026-10-21: boss 22 -> 29 mm so that the handle turns in front of the dial (57.0) and the hands (62.9):
# seated on the drum shoulder (z 41.6) the handle now spans z 64.6..70.6. Pocket depth (13 mm) unchanged.
KEY_L = 29.0
key = cq.Solid.makeBox(60, 12, 6, V(KX - 30, KY - 6, 0)).fuse(cyl(KX, KY, KEY_R, 0, KEY_L)).cut(prism_poly(hexagon(HEX_KEY), KEY_L - 13, KEY_L + 1, at=(KX, KY)))
add("winding_key", "巻きカギ（六角10 mm）", key, group="accessory", color="accent")

# ------------------------------------------------------------------ outputs
def print_oriented(p):
    s = p["solid"]
    if p["print_rot"]:
        ax, ang = p["print_rot"]
        s = s.rotate(V(0, 0, 0), V(1, 0, 0) if ax == "x" else V(0, 1, 0), ang)
    bb = s.BoundingBox()
    return s.translate(V(-(bb.xmin + bb.xmax) / 2, -(bb.ymin + bb.ymax) / 2, -bb.zmin))

COLORS = {"plate": (0.84, 0.86, 0.84, 1), "pillar": (0.55, 0.57, 0.58, 1), "steel": (0.72, 0.74, 0.78, 1), "gear1": (0.80, 0.58, 0.20, 1),
          "gear2": (0.86, 0.66, 0.26, 1), "gear3": (0.90, 0.74, 0.36, 1), "gear4": (0.78, 0.52, 0.22, 1), "escape": (0.72, 0.30, 0.22, 1),
          "anchor": (0.20, 0.36, 0.42, 1), "accent": (0.22, 0.45, 0.50, 1), "drum": (0.45, 0.33, 0.24, 1), "hand": (0.10, 0.11, 0.13, 1),
          "dial": (0.95, 0.94, 0.90, 1), "hammer": (0.20, 0.36, 0.42, 1), "bell": (0.85, 0.72, 0.35, 1), "bob": (0.80, 0.58, 0.20, 1),
          "weight": (0.33, 0.35, 0.38, 1)}

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    print(len(PARTS), "components")
    for p in PARTS:
        s = p["solid"]
        assert s.isValid(), p["name"]
        p["vol"] = s.Volume()
        bb = s.BoundingBox()
        p["bb"] = (bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax)
    if mode in ("all", "check"):
        bad = 0
        skip = lambda p: p["name"].startswith(("rod_", "tie_rod", "weight_pin", "pendulum_rod"))
        for a, b in itertools.combinations(PARTS, 2):
            if skip(a) or skip(b): continue
            ba, bb_ = a["bb"], b["bb"]
            if any(ba[i] > bb_[i + 3] or bb_[i] > ba[i + 3] for i in range(3)): continue
            try:
                v = a["solid"].intersect(b["solid"]).Volume()
            except Exception as ex:
                v = -1
            if abs(v) > 0.3:
                bad += 1
                print(f"  INTERFERENCE {a['name']} x {b['name']}: {v:.2f} mm3")
        print("interference pairs:", bad)
    if mode in ("all", "export"):
        asm = cq.Assembly(name="pendulum_clock")
        seen = set()
        bom = []
        for p in PARTS:
            asm.add(p["solid"], name=p["name"], color=cq.Color(*COLORS[p["color"]]))
            key_ = p["share"] or p["name"]
            if p["kind"] == "print" and key_ not in seen:
                seen.add(key_)
                s = print_oriented(p)
                cq.exporters.export(cq.Workplane(obj=s), f"{OUT}/parts_step/{key_}.step")
                cq.exporters.export(cq.Workplane(obj=s), f"{OUT}/parts_stl/{key_}.stl", tolerance=0.02, angularTolerance=0.15)
                bb = s.BoundingBox()
                bom.append(dict(file=key_, ja=p["ja"], size=[round(bb.xlen, 1), round(bb.ylen, 1), round(bb.zlen, 1)], vol=round(p["vol"] / 1000, 1)))
        os.makedirs(f"{OUT}/extras", exist_ok=True)
        for nm, sh in ANCHOR_VARIANTS.items():
            DV = D.variant(A_DIST="RW * sqrt(2) + %g" % sh)
            sv = print_oriented(dict(solid=anchor_solid(DV.anchor_polygon()), print_rot=None))
            cq.exporters.export(cq.Workplane(obj=sv), f"{OUT}/extras/anchor_{nm}.step")
            cq.exporters.export(cq.Workplane(obj=sv), f"{OUT}/extras/anchor_{nm}.stl", tolerance=0.02, angularTolerance=0.15)
        asm.save(f"{OUT}/pendulum_clock_assembly.step")
        json.dump(bom, open(f"{OUT}/bom_print.json", "w"), ensure_ascii=False, indent=1)
        print("exported", len(bom), "printable part files")
    if mode in ("all", "mesh", "export"):
        blob = bytearray(); meta = []; ntri = 0
        for p in PARTS:
            if p["group"] == "accessory": continue
            small = p["vol"] < 600
            vs, ts = p["solid"].tessellate(0.12 if small else 0.06, 0.5 if small else 0.28)
            v = np.array([[q.x, q.y, q.z] for q in vs], dtype=np.float64)
            t = np.array(ts, dtype=np.uint32)
            lo, hi = v.min(0), v.max(0)
            sc = np.maximum(hi - lo, 1e-6)
            q = np.round((v - lo) / sc * 65535).astype(np.uint16)
            off = len(blob); blob += q.tobytes()
            idx16 = len(v) < 65536
            ioff = len(blob); blob += (t.astype(np.uint16) if idx16 else t).tobytes()
            if len(blob) % 4: blob += b"\0" * (4 - len(blob) % 4)
            ntri += len(t)
            meta.append(dict(n=p["name"], ja=p["ja"], g=p["group"], c=p["color"], k=p["kind"], lo=[round(float(x), 4) for x in lo],
                             sc=[round(float(x), 4) for x in sc], nv=len(v), nt=len(t), o=off, io=ioff, i16=idx16))
        open(f"{OUT}/mesh.bin", "wb").write(blob)
        json.dump(meta, open(f"{OUT}/mesh.json", "w"), ensure_ascii=False)
        print("mesh: %d tris, %.2f MB" % (ntri, len(blob) / 1e6))
