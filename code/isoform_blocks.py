from difflib import SequenceMatcher

MINBLOCK = 5

def blocks(ref, iso, minblock=MINBLOCK):
    sm = SequenceMatcher(None, ref, iso, autojunk=False)
    return [(b.a, b.b, b.size) for b in sm.get_matching_blocks() if b.size >= minblock]

def segments(blocks_, ref_len, iso_len):
    nov, dele = [], []
    pa, pb = 0, 0
    for k, (a, b, s) in enumerate(blocks_):
        where = "N" if k == 0 else "I"
        if b > pb: nov.append((pb, b, where))
        if a > pa: dele.append((pa, a, where))
        pa, pb = a + s, b + s
    if iso_len > pb: nov.append((pb, iso_len, "C" if blocks_ else "N"))
    if ref_len > pa: dele.append((pa, ref_len, "C" if blocks_ else "N"))
    return nov, dele

def classify(ref, iso, minblock=MINBLOCK):
    bl = blocks(ref, iso, minblock)
    nov, dele = segments(bl, len(ref), len(iso))
    shared = sum(s for _, _, s in bl)
    n_novel, n_del = len(iso) - shared, len(ref) - shared
    novN = sum(e - s for s, e, w in nov if w == "N"); novC = sum(e - s for s, e, w in nov if w == "C")
    novI = sum(e - s for s, e, w in nov if w == "I")
    delN = sum(e - s for s, e, w in dele if w == "N"); delC = sum(e - s for s, e, w in dele if w == "C")
    delI = sum(e - s for s, e, w in dele if w == "I")
    if not bl:
        cls, sub = "unrelated", "no shared block"
    elif n_novel == 0 and n_del == 0:
        cls, sub = "identical", ""
    elif n_novel == 0:
        if delI == 0:
            cls = "truncation"
            sub = "N-terminal" if delN and not delC else ("C-terminal" if delC and not delN else "both termini")
        else:
            cls = "internal deletion"
            sub = "internal only" if not (delN or delC) else "internal plus terminal"
    elif n_del == 0:
        cls = "addition"
        parts = [p for p, v in (("N-terminal extension", novN), ("C-terminal extension", novC), ("internal insertion", novI)) if v]
        sub = " + ".join(parts)
    else:
        cls = "replacement"
        parts = []
        if novN or delN: parts.append("N-terminal")
        if novC or delC: parts.append("C-terminal")
        if novI or delI: parts.append("internal")
        sub = " + ".join(parts)
    rel = "longer" if len(iso) > len(ref) else ("equal" if len(iso) == len(ref) else "shorter")
    longest = max([s for _, _, s in bl], default=0)
    return {
        "class": cls, "subclass": sub, "length_relation": rel,
        "ref_len": len(ref), "iso_len": len(iso), "shared": shared,
        "novel": n_novel, "deleted": n_del,
        "novel_N": novN, "novel_C": novC, "novel_I": novI,
        "deleted_N": delN, "deleted_C": delC, "deleted_I": delI,
        "n_blocks": len(bl), "longest_block": longest,
        "novel_segments": ";".join(f"{s+1}-{e}:{w}" for s, e, w in nov),
        "deleted_segments": ";".join(f"{s+1}-{e}:{w}" for s, e, w in dele),
    }
