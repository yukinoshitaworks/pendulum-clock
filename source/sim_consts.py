"""Write the design values the simulator needs into its page (the `const S = {...}` line).

    python3 sim_consts.py ../index.html        # updates the file in place

Only the keys that come from design.py are replaced; the rest of S (strike table, positions, etc.) is kept.
"""
import sys, json, re
import design as D

def values():
    r5 = lambda v: round(float(v), 5)
    return {
        "tooth": [[r5(x), r5(y)] for x, y in D.esc_tooth_sim()],
        "pallets": [[[r5(x), r5(y)] for x, y in pc] for pc in D.pallet_pieces()],
        "escPhase": float(D.ESC_PHASE), "pitch": float(D.PITCH), "nEsc": D.N_ESC, "aDist": float(D.A_DIST),
        "lift": float(D.LIFT), "lock": float(D.LOCK), "drop": float(D.DROP),
        "rThird": D.R_THIRD, "rCentre": D.R_CENTRE, "rGreat": D.R_GREAT, "rMinw": D.R_MINW, "rHour": D.R_HOUR,
        "cordR": D.CORD_R, "drumTurns": D.DRUM_TURNS, "weightX": D.WEIGHT_X, "weightZ": D.WEIGHT_Z,
        "teeth": {"great": D.great.z, "cPin": D.c_pin.z, "cWhl": D.c_whl.z, "tPin": D.t_pin.z, "tWhl": D.t_whl.z,
                  "ePin": D.e_pin.z, "cannon": D.cannon.z, "mWhl": D.m_whl.z, "mPin": D.m_pin.z, "hWhl": D.h_whl.z},
    }

if __name__ == "__main__":
    path = sys.argv[1]
    lines = open(path, encoding="utf-8").read().split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith("const S = {"))
    S = json.loads(lines[i][len("const S = "):].rstrip(";"))
    S.update(values())
    lines[i] = "const S = " + json.dumps(S, ensure_ascii=False) + ";"
    open(path, "w", encoding="utf-8").write("\n".join(lines))
    print("updated", path, "rGreat 1/%.0f cordR %.1f lift %.2f lock %.2f drop %.2f" %
          (1 / S["rGreat"], S["cordR"], S["lift"] * 57.2958, S["lock"] * 57.2958, S["drop"] * 57.2958))
