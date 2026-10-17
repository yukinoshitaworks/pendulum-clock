"""Weight-driven pendulum clock -- single source of truth for CAD and simulation.
Units: mm, radians. World frame: x right, y up, z toward the viewer (front).
z = 0 is the inner face of the back plate.
"""
import math
from math import pi, sin, cos, tan, atan2, acos, sqrt, radians, degrees
import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity
from shapely.ops import unary_union

ALPHA = radians(20.0)
BACKLASH = 0.12          # tooth thinning per gear at the pitch circle
PLATE_T = 5.0
GAP = 36.0               # distance between plate inner faces
ROD = 3.0                # arbor rod diameter

# ---------------------------------------------------------------- layout
C = (0.0, 0.0)                       # centre arbor (minute hand)
G = (0.0, -63.0)                     # great wheel / drum
E = (0.0, 60.0)                      # escape wheel
_yT = (52**2 - 44**2 + 60**2) / (2 * 60)
T = (-sqrt(52**2 - _yT**2), _yT)     # third arbor
RW = 30.0                            # escape wheel tip radius
N_ESC = 30
A_DIST = RW * sqrt(2)                # 7.5 tooth span
A = (0.0, E[1] + A_DIST)             # anchor / pendulum pivot
M = (0.0, -24.0)                     # minute wheel stud
H = (-T[0], T[1])                    # hammer stud
P = [(-38.0, 116.0), (38.0, 116.0), (-58.0, -30.0), (58.0, -30.0)]  # pillars
D0 = (0.0, 79.0)                     # dial top standoff

Z_FRONT_IN = GAP
Z_FRONT_OUT = GAP + PLATE_T          # 41

# ---------------------------------------------------------------- gears
def inv(a):
    return tan(a) - a


class Gear:
    def __init__(self, name, z, m, x, center, z0, z1):
        self.name, self.z, self.m, self.x = name, z, m, x
        self.center, self.z0, self.z1 = center, z0, z1
        self.r = m * z / 2
        self.rb = self.r * cos(ALPHA)
        self.ra = self.r + m * (1 + x)
        self.rf = self.r - m * (1.25 - x)
        self.s = m * (pi / 2 + 2 * x * tan(ALPHA)) - BACKLASH
        self.phase = 0.0

    def half(self, rho):
        rho = max(rho, self.rb)
        return self.s / (2 * self.r) + inv(ALPHA) - inv(acos(self.rb / rho))

    def limit_tip(self, min_land):
        while self.half(self.ra) * 2 * self.ra < min_land:
            self.ra -= 0.01

    def flank(self, n):
        """(rho, half_angle) samples from flank start to tip."""
        rs = max(self.rb, self.rf)
        t0 = sqrt(max(rs**2 - self.rb**2, 0)) / self.rb
        t1 = sqrt(self.ra**2 - self.rb**2) / self.rb
        out = []
        for i in range(n):
            t = t0 + (t1 - t0) * i / (n - 1)
            rho = self.rb * sqrt(1 + t * t)
            out.append((rho, self.half(rho)))
        return out

    def segments(self, n=6):
        """Closed outline as CAD-friendly segments, gear-local coords, tooth 0 at self.phase."""
        segs = []
        fl = self.flank(n)
        pol = lambda rho, a: (rho * cos(a), rho * sin(a))
        pitch = 2 * pi / self.z
        radial = self.rf < self.rb - 1e-9
        h_root = fl[0][1]
        for k in range(self.z):
            th = self.phase + k * pitch
            if radial:
                segs.append(("line", [pol(self.rf, th - h_root), pol(fl[0][0], th - h_root)]))
            segs.append(("spline", [pol(rho, th - h) for rho, h in fl]))
            ht = fl[-1][1]
            segs.append(("arc", [pol(self.ra, th - ht), pol(self.ra, th), pol(self.ra, th + ht)]))
            segs.append(("spline", [pol(rho, th + h) for rho, h in reversed(fl)]))
            if radial:
                segs.append(("line", [pol(fl[0][0], th + h_root), pol(self.rf, th + h_root)]))
            a0, a1 = th + h_root, th + pitch - h_root
            segs.append(("arc", [pol(self.rf, a0), pol(self.rf, (a0 + a1) / 2), pol(self.rf, a1)]))
        return segs

    def polygon(self, n=10):
        pts = []
        for kind, p in self.segments(n):
            if kind == "arc":
                a0 = atan2(p[0][1], p[0][0]); a1 = atan2(p[2][1], p[2][0])
                while a1 < a0: a1 += 2 * pi
                rr = math.hypot(*p[0])
                m_ = max(2, int((a1 - a0) / radians(3)) + 1)
                pts += [(rr * cos(a0 + (a1 - a0) * i / m_), rr * sin(a0 + (a1 - a0) * i / m_)) for i in range(m_)]
            else:
                pts += p[:-1]
        return Polygon(pts)


def pair(wheel, pinion):
    """Tune tips/roots of a meshing pair and set phases so teeth interleave at the reference pose."""
    a = math.dist(wheel.center, pinion.center)
    pinion.limit_tip(0.35 * pinion.m)
    lim = sqrt(wheel.rb**2 + (a * sin(ALPHA))**2) - 0.03       # avoid interference below pinion base circle
    wheel.ra = min(wheel.ra, lim)
    wheel.rf = a - pinion.ra - 0.25 * wheel.m
    pinion.rf = a - wheel.ra - 0.25 * wheel.m
    beta = atan2(pinion.center[1] - wheel.center[1], pinion.center[0] - wheel.center[0])
    wheel.phase = beta                                 # a tooth points at the pinion
    pinion.phase = beta + pi + pi / pinion.z           # a gap points back
    eps = (sqrt(pinion.ra**2 - pinion.rb**2) + sqrt(wheel.ra**2 - wheel.rb**2) - a * sin(ALPHA)) / (pi * wheel.m * cos(ALPHA))
    return a, eps


great = Gear("great", 72, 1.5, -0.30, G, 9.0, 14.0)
c_pin = Gear("centre pinion", 12, 1.5, 0.30, C, 7.0, 15.0)
c_whl = Gear("centre wheel", 96, 1.0, -0.45, C, 2.0, 7.0)
t_pin = Gear("third pinion", 8, 1.0, 0.45, T, 1.5, 15.5)
t_whl = Gear("third wheel", 80, 1.0, -0.45, T, 15.5, 20.5)
e_pin = Gear("escape pinion", 8, 1.0, 0.45, E, 14.0, 22.0)
cannon = Gear("cannon pinion", 12, 1.0, 0.30, C, 41.6, 46.6)
m_whl = Gear("minute wheel", 36, 1.0, -0.30, M, 42.0, 46.0)
m_pin = Gear("minute pinion", 10, 0.96, 0.40, M, 46.0, 53.0)
h_whl = Gear("hour wheel", 40, 0.96, -0.40, C, 47.5, 52.0)

PAIRS = [(great, c_pin), (c_whl, t_pin), (t_whl, e_pin), (m_whl, cannon), (h_whl, m_pin)]
PAIR_INFO = [pair(w, p) for w, p in PAIRS]

# train ratios relative to escape-wheel advance a (one tooth pitch = one pendulum period)
R_THIRD = e_pin.z / t_whl.z                 # 1/10
R_CENTRE = R_THIRD * t_pin.z / c_whl.z      # 1/120
R_GREAT = R_CENTRE * c_pin.z / great.z      # 1/720
R_MINW = R_CENTRE * cannon.z / m_whl.z
R_HOUR = R_MINW * m_pin.z / h_whl.z         # 1/1440


def check_pair(wheel, pinion, steps=48):
    """Rotate the pair through two pitches; return (max overlap area, min gap)."""
    pw, pp = wheel.polygon(), pinion.polygon()
    worst, gap = 0.0, 1e9
    for i in range(steps):
        th = 2 * (2 * pi / wheel.z) * i / steps
        a = affinity.translate(affinity.rotate(pw, th, origin=(0, 0), use_radians=True), *wheel.center)
        b = affinity.translate(affinity.rotate(pp, -th * wheel.z / pinion.z, origin=(0, 0), use_radians=True), *pinion.center)
        worst = max(worst, a.intersection(b).area)
        gap = min(gap, a.distance(b))
    return worst, gap


# ---------------------------------------------------------------- escapement (Graham deadbeat)
LIFT = radians(3.0)       # anchor rotation during impulse on each pallet
LOCK = radians(1.5)
DROP = radians(2.0)       # escape-wheel free travel between pallets
PITCH = 2 * pi / N_ESC
PH_UE = (LIFT - LOCK) / 2       # entry unlock (anchor angle, CCW +)
PH_LE = -(LIFT + LOCK) / 2      # entry let-off
PH_UX = (LOCK - LIFT) / 2       # exit unlock
PH_LX = (LIFT + LOCK) / 2       # exit let-off
E_LOC = (0.0, -A_DIST)          # wheel centre in anchor frame


def rot(p, a):
    return (p[0] * cos(a) - p[1] * sin(a), p[0] * sin(a) + p[1] * cos(a))


def _isect(rho, side):
    y = (RW**2 - A_DIST**2 - rho**2) / (2 * A_DIST)
    return (side * sqrt(rho**2 - y**2), y)


def _travel(p):  # CW travel coordinate of a point about the wheel centre
    return -atan2(p[1] - E_LOC[1], p[0] - E_LOC[0])


def _solve_t():
    lo, hi = 0.5, 4.0
    for _ in range(60):
        t = (lo + hi) / 2
        d = (_travel(_isect(RW - t / 2, 1)) - _travel(_isect(RW - t / 2, -1))) % PITCH
        if d > DROP: lo = t
        else: hi = t
    return t


PALLET_T = _solve_t()
RI, RO = RW - PALLET_T / 2, RW + PALLET_T / 2
P_EO, P_EI = _isect(RO, -1), _isect(RI, -1)
P_XI, P_XO = _isect(RI, 1), _isect(RO, 1)
ENT_LOCK = rot(P_EO, -PH_UE); ENT_LET = rot(P_EI, -PH_LE)
EXI_LOCK = rot(P_XI, -PH_UX); EXI_LET = rot(P_XO, -PH_LX)
NIB_SPAN = radians(17)


def _arc(rho, a0, a1, n):
    return [(rho * cos(a0 + (a1 - a0) * i / n), rho * sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]


def pallet_pieces(n=4):
    """Convex slices of both pallet nibs in the anchor frame (for collision in the simulation)."""
    out = []
    for lock, let, sgn in ((ENT_LOCK, ENT_LET, -1), (EXI_LOCK, EXI_LET, 1)):
        outer, inner = (lock, let) if sgn < 0 else (let, lock)
        ao, ai = atan2(outer[1], outer[0]), atan2(inner[1], inner[0])
        back = (radians(225) - NIB_SPAN) if sgn < 0 else (radians(315) + NIB_SPAN)
        ao %= 2 * pi; ai %= 2 * pi
        o = _arc(RO, ao, back, n); i = _arc(RI, ai, back, n)
        for k in range(n):
            out.append([o[k], o[k + 1], i[k + 1], i[k]])
    return out


def anchor_polygon():
    parts = [Point(0, 0).buffer(6.5, 12)]
    for sgn in (-1, 1):
        back = (radians(225) - NIB_SPAN) if sgn < 0 else (radians(315) + NIB_SPAN)
        lock, let = (ENT_LOCK, ENT_LET) if sgn < 0 else (EXI_LOCK, EXI_LET)
        outer, inner = (lock, let) if sgn < 0 else (let, lock)
        ao, ai = atan2(outer[1], outer[0]) % (2 * pi), atan2(inner[1], inner[0]) % (2 * pi)
        back2 = back - sgn * radians(4)
        nib = Polygon(_arc(RO, ao, back2, 12) + list(reversed(_arc(RI, ai, back2, 12))))
        end = (RW * cos(back), RW * sin(back))
        arm = LineString([(0, 0), (end[0] * 0.5, end[1] * 0.5 + 3.0), end]).buffer(3.0, 6)
        parts += [nib, arm, Point(end).buffer(3.6, 8)]
    return unary_union(parts)


# escape tooth (wheel-local, leading corner of tooth 0 at angle 0 on the tip circle; wheel turns CW)
ESC_ROOT = 24.0
ESC_LAND = 0.5 / RW
ESC_UNDERCUT = radians(8)
ESC_BACK = radians(20)


def esc_tooth():
    h = RW - ESC_ROOT
    front_root = h * tan(ESC_UNDERCUT) / ESC_ROOT
    back_root = ESC_LAND + h * tan(ESC_BACK) / ESC_ROOT
    pol = lambda rho, a: (rho * cos(a), rho * sin(a))
    return [pol(RW, 0), pol(RW, ESC_LAND), pol(ESC_ROOT, back_root), pol(ESC_ROOT, front_root)]


# phase: at the reference pose (a = 0, anchor centred) one tooth rests on the entry impulse face
def _tip_on_entry_face():
    dx, dy = ENT_LET[0] - ENT_LOCK[0], ENT_LET[1] - ENT_LOCK[1]
    fx, fy = ENT_LOCK[0] - E_LOC[0], ENT_LOCK[1] - E_LOC[1]
    qa, qb, qc = dx * dx + dy * dy, 2 * (fx * dx + fy * dy), fx * fx + fy * fy - RW**2
    t = [(-qb + s_ * sqrt(qb * qb - 4 * qa * qc)) / (2 * qa) for s_ in (-1, 1)]
    t = [v for v in t if -0.01 <= v <= 1.01][0]
    return atan2(fy + t * dy, fx + t * dx)


ESC_PHASE = _tip_on_entry_face() + radians(0.05)


def esc_polygon():
    t = esc_tooth(); pts = []
    for k in range(N_ESC):
        a = ESC_PHASE + k * PITCH
        pts += [rot(t[3], a), rot(t[0], a), rot(t[1], a), rot(t[2], a)]
    return Polygon(pts)


# ---------------------------------------------------------------- strike: snail cam + gravity hammer
CAM_RMIN, CAM_RMAX = 6.0, 16.0
CAM_Z = (24.0, 30.0)
_ch = math.dist(C, H)
_bc = atan2(H[1] - C[1], H[0] - C[0]) - acos(11.0 / _ch)        # contact direction where follower motion is radial
F0 = (C[0] + 11 * cos(_bc), C[1] + 11 * sin(_bc))
TAIL_U = ((F0[0] - H[0]) / math.dist(F0, H), (F0[1] - H[1]) / math.dist(F0, H))   # tail direction at gamma = 0
TAIL_S0, TAIL_S1, TAIL_W = 4.0, sqrt(_ch**2 - CAM_RMAX**2) - 1.0, 5.0
TAIL_N = (TAIL_U[1], -TAIL_U[0])        # unit normal pointing toward the cam side... fixed below
if (C[0] - H[0]) * TAIL_N[0] + (C[1] - H[1]) * TAIL_N[1] < 0:
    TAIL_N = (-TAIL_N[0], -TAIL_N[1])
HEAD_L, HEAD_R, BELL_R, BELL_GAP = 50.0, 7.0, 30.0, 0.8
TAIL_RAKE = tan(radians(35))      # end face raked back so it clears the cam step while falling


def cam_r(alpha):
    return CAM_RMIN + (CAM_RMAX - CAM_RMIN) * (alpha % (2 * pi)) / (2 * pi)


def tail_samples():
    """Tail contact edge + end face sample points, hammer-local (origin H, gamma = 0)."""
    pts = []
    s = TAIL_S0
    while s < TAIL_S1:
        pts.append((TAIL_U[0] * s, TAIL_U[1] * s)); s += 0.25
    for w in np.linspace(0, TAIL_W, 9):
        sw = TAIL_S1 - w * TAIL_RAKE
        pts.append((TAIL_U[0] * sw - TAIL_N[0] * w, TAIL_U[1] * sw - TAIL_N[1] * w))
    return pts


_TS = np.array(tail_samples())


def hammer_gamma(cam_rot):
    """Hammer angle (CCW +, lift direction) resting on the cam whose local frame is rotated by cam_rot."""
    def pen(g):
        c, s = cos(g), sin(g)
        x = H[0] + _TS[:, 0] * c - _TS[:, 1] * s - C[0]
        y = H[1] + _TS[:, 0] * s + _TS[:, 1] * c - C[1]
        al = (np.arctan2(y, x) - cam_rot) % (2 * pi)
        rr = CAM_RMIN + (CAM_RMAX - CAM_RMIN) * al / (2 * pi)
        return bool(np.any(np.hypot(x, y) < rr))
    lo, hi = radians(-12), radians(12)
    for _ in range(30):
        mid = (lo + hi) / 2
        if pen(mid): lo = mid
        else: hi = mid
    return hi


def build_hammer_table(n=1440):
    """gamma vs centre-arbor angle theta_c (CW +). Cam orientation chosen so the drop happens at theta_c = 0."""
    raw = [hammer_gamma(-2 * pi * i / n) for i in range(n)]
    d = [raw[(i + 1) % n] - raw[i] for i in range(n)]
    k = int(np.argmin(d)) + 1                 # first sample after the drop
    cam0 = -2 * pi * k / n                    # rotate cam so the drop sits just before theta_c = 0
    tab = [hammer_gamma(cam0 - 2 * pi * i / n) for i in range(n)]
    return cam0, tab


CAM0, HAMMER_TAB = build_hammer_table()
G_REST, G_MAX = HAMMER_TAB[0], max(HAMMER_TAB)
_ha = radians(14.0) + G_REST
HEAD_C = (H[0] + HEAD_L * cos(_ha), H[1] + HEAD_L * sin(_ha))      # head centre at rest
_v = (sin(_ha), -cos(_ha))                                         # head velocity direction when falling
BELL = (HEAD_C[0] + (HEAD_R + BELL_GAP + BELL_R) * _v[0], HEAD_C[1] + (HEAD_R + BELL_GAP + BELL_R) * _v[1])


def cam_polygon(n=180):
    """Cam outline in world orientation at theta_c = 0 (centre C at origin)."""
    pts = [(cam_r(2 * pi * i / n) * cos(CAM0 + 2 * pi * i / n), cam_r(2 * pi * i / n) * sin(CAM0 + 2 * pi * i / n)) for i in range(n)]
    pts.append((CAM_RMAX * cos(CAM0), CAM_RMAX * sin(CAM0)))
    return Polygon(pts)


def hammer_polygon():
    """Hammer outline at rest pose, world coords relative to H."""
    g = G_REST
    u, nn = rot(TAIL_U, g), rot(TAIL_N, g)
    tail = Polygon([(u[0] * 0, u[1] * 0), (u[0] * TAIL_S1, u[1] * TAIL_S1),
                    (u[0] * (TAIL_S1 - TAIL_W * TAIL_RAKE) - nn[0] * TAIL_W, u[1] * (TAIL_S1 - TAIL_W * TAIL_RAKE) - nn[1] * TAIL_W),
                    (-nn[0] * TAIL_W, -nn[1] * TAIL_W)])
    hc = (HEAD_C[0] - H[0], HEAD_C[1] - H[1])
    arm = LineString([(0, 0), hc]).buffer(1.9, 4)
    return unary_union([tail, arm, Point(0, 0).buffer(6.0, 12), Point(hc).buffer(HEAD_R, 12)])


# ---------------------------------------------------------------- pendulum / weight
PEND_T = 1.0
PEND_L = 9806.65 * (PEND_T / (2 * pi))**2     # 248.4 mm equivalent length
BOB_Y = 250.0                                 # bob centre below pivot (rating nut adjusts)
PEND_Z = -13.75
DRUM_R = 20.0
CORD_R = DRUM_R + 0.5
DRUM_TURNS = 9.0
WEIGHT_X = G[0] - CORD_R
WEIGHT_Z = 26.5

if __name__ == "__main__":
    print("T", T, "A", A, "H", H)
    for (w, p), (a, eps) in zip(PAIRS, PAIR_INFO):
        ov, gap = check_pair(w, p)
        print(f"{w.name:13s} {w.z:3d} / {p.name:14s} {p.z:2d}  a={a:6.3f} eps={eps:4.2f} "
              f"ra_w={w.ra:6.2f} rf_w={w.rf:6.2f} ra_p={p.ra:5.2f} rf_p={p.rf:5.2f} rb_p={p.rb:5.2f} overlap={ov:.4f} mingap={gap:.3f}")
    print("ratios", 1 / R_THIRD, 1 / R_CENTRE, 1 / R_GREAT, 1 / R_HOUR)
    print("pallet t", PALLET_T, "RI", RI, "RO", RO)
    print("cam0", degrees(CAM0), "rest", degrees(G_REST), "max", degrees(G_MAX), "lift", degrees(G_MAX - G_REST))
    print("head", HEAD_C, "bell", BELL, "tail_s1", TAIL_S1)
    d = np.diff(HAMMER_TAB)
    print("monotone rise except drop:", int((d < -1e-4).sum()), "neg steps; biggest fall", degrees(d.min()), "at", int(d.argmin()))
