"""Energy budget of the escapement and pendulum, with friction (the simulator has none).

    python3 escapement_budget.py

Per pendulum period (1 s) the escapement gives two impulses and has the tooth resting on a locking face for the
rest of each beat. In a deadbeat the locking faces are arcs about the anchor pivot, so while a tooth rests on one
the pendulum drags it through the supplementary arc against friction  mu * F * r_lock.  Both the impulse and that
drag are proportional to the drive torque, so their balance gives an amplitude that more weight cannot raise:

    theta_cap = LOCK/2 + dphi_imp * eta_imp / (2 * mu * r_lock / RW)

dphi_imp  escape-wheel turn while a tooth is on an impulse face (= half a pitch - drop)
eta_imp   impulse-face efficiency, a wedge: tan(a) / tan(a + atan(mu)), tan(a) = LIFT * r_lock / (RW * dphi_imp)
"""
import math
from math import radians, degrees, atan, tan, sqrt, pi
import design as D

G0 = 9.80665


def budget(mu, mass=1.5, eta_train=0.6, drop=D.DROP, lift=D.LIFT, lock=D.LOCK, r_lock_ratio=None,
           pend_I=None, Q_air=150.0, T_bearing=0.04):
    """Return dict with cap amplitude, running amplitude (pendulum losses included) and margins, all in degrees.
    T_bearing: Coulomb friction torque of the two anchor bearings under the pendulum load (N*mm)."""
    RW = D.RW
    r_lock = sqrt(D.A_DIST**2 - RW**2) if r_lock_ratio is None else r_lock_ratio * RW
    dphi = D.PITCH / 2 - drop
    ta = lift * r_lock / (RW * dphi)
    a = atan(ta)
    eta_imp = ta / tan(a + atan(mu))
    mu_l = mu * r_lock / RW
    letoff = (lift + lock) / 2
    cap = lock / 2 + dphi * eta_imp / (2 * mu_l)
    T_esc = mass * G0 * D.CORD_R * D.R_GREAT * eta_train          # N*mm at the escape arbor
    if pend_I is None:
        pend_I = pendulum_inertia()
    w = 2 * pi / D.PEND_T

    def e_in(th):   # J per period
        return 2 * T_esc * 1e-3 * (dphi * eta_imp - mu_l * (2 * th - lock))

    def e_loss(th):
        return 2 * pi * 0.5 * pend_I * w * w * th * th / Q_air + 4 * T_bearing * 1e-3 * th

    # running amplitude: largest th > letoff with e_in >= e_loss
    run = None
    th = cap
    while th > letoff:
        if e_in(th) >= e_loss(th):
            run = th
            break
        th -= radians(0.005)
    return dict(mu=mu, eta_imp=eta_imp, letoff=degrees(letoff), cap=degrees(cap),
                run=None if run is None else degrees(run),
                margin=None if run is None else degrees(run - letoff),
                e_in_per_period_uJ=e_in(letoff + radians(0.5)) * 1e6, T_esc=T_esc)


def pendulum_inertia(bob_g=80.0, rod_g=22.0):
    """kg*m^2 about the pivot: bob at BOB_Y, M4 rod from 15 to 300 mm below the pivot."""
    L = D.BOB_Y / 1000
    l1, l2 = 0.015, 0.300
    return bob_g / 1000 * L * L + rod_g / 1000 * (l2**3 - l1**3) / (3 * (l2 - l1))


def row(name, **kw):
    out = [name]
    for mu in (0.15, 0.25, 0.35, 0.45):
        b = budget(mu, **kw)
        run = "stops" if b["run"] is None else "%.2f (+%.2f)" % (b["run"], b["margin"])
        out.append("%5.2f / %-13s" % (b["cap"], run))
    print("%-34s" % out[0] + "  ".join(out[1:]))


if __name__ == "__main__":
    b = budget(0.35)
    print("escape-arbor torque (1.5 kg, train eff 0.6): %.2f N*mm; let-off amplitude %.2f deg" % (b["T_esc"], b["letoff"]))
    print("pendulum inertia %.2e kg*m^2 (bob 80 g + rod 22 g)" % pendulum_inertia())
    print()
    print("amplitude in degrees: cap (friction only, any weight) / running with pendulum losses (+ margin over let-off)")
    print("%-34s" % "" + "  ".join("mu = %-16.2f" % m for m in (0.15, 0.25, 0.35, 0.45)))
    row("current 1.5 kg")
    row("current 3.0 kg", mass=3.0)
    row("drop 2.75 -> 2.0", drop=radians(2.0))
    row("drop 2.75 -> 1.5", drop=radians(1.5))
    row("lift 1.5 -> 2.5", lift=radians(2.5))
    row("span 7.5 -> 5.5 teeth (r_lock .65)", r_lock_ratio=tan(radians(5.5 * 12 / 2)))
    row("pendulum losses x2", Q_air=75.0, T_bearing=0.08)
