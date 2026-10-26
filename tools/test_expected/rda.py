"""Minimal reader of R .rda files (XDR serialization, versions 2 and 3).

Supports the object types found in data-frame datasets: vectors (logical,
integer, double, character), lists, pairlists, symbols, attributes, and the
compact ALTREP sequences used for row names.
"""
import gzip, bz2, lzma, struct
import numpy as np

NA_INT = -2147483648


class RObj:
    def __init__(self, typ, value, attr=None):
        self.typ, self.value, self.attr = typ, value, attr or {}

    def __repr__(self):
        return f"RObj({self.typ}, attr={list(self.attr)})"


class Reader:
    def __init__(self, buf):
        self.b, self.p, self.refs = buf, 0, []

    def int(self):
        v = struct.unpack(">i", self.b[self.p:self.p + 4])[0]
        self.p += 4
        return v

    def length(self):
        n = self.int()
        if n == -1:
            hi, lo = self.int(), self.int()
            n = (hi << 32) + lo
        return n

    def attrs(self, has_attr):
        if not has_attr:
            return {}
        pl = self.item()
        return dict(pl.value) if pl is not None else {}

    def item(self):
        flags = self.int()
        typ = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))
        if typ == 254:  # NILVALUE
            return None
        if typ in (253, 252, 251, 250, 242, 241):
            return RObj("special", typ)
        if typ == 255:  # REFSXP
            idx = flags >> 8
            if idx == 0:
                idx = self.int()
            return self.refs[idx - 1]
        if typ == 1:  # SYMSXP
            s = self.item()
            sym = RObj("sym", s.value)
            self.refs.append(sym)
            return sym
        if typ in (2, 6, 239, 240):  # pairlist / language
            out = []
            node_attr = self.attrs(has_attr)
            while True:
                tag = self.item().value if has_tag else None
                car = self.item()
                out.append((tag, car))
                flags = self.int()
                typ2 = flags & 0xFF
                if typ2 == 254:
                    break
                if typ2 not in (2, 6, 239, 240):
                    raise ValueError(f"unexpected cdr type {typ2}")
                has_attr = bool(flags & (1 << 9))
                has_tag = bool(flags & (1 << 10))
                if has_attr:
                    self.item()
            return RObj("pairlist", out, node_attr)
        if typ == 9:  # CHARSXP
            n = self.int()
            if n == -1:
                return RObj("char", None)
            s = self.b[self.p:self.p + n].decode("utf-8", "replace")
            self.p += n
            return RObj("char", s)
        if typ in (10, 13):  # LGLSXP, INTSXP
            n = self.length()
            v = np.frombuffer(self.b, ">i4", n, self.p).astype(np.int64)
            self.p += 4 * n
            return RObj("int" if typ == 13 else "lgl", v, self.attrs(has_attr))
        if typ == 14:  # REALSXP
            n = self.length()
            v = np.frombuffer(self.b, ">f8", n, self.p).astype(float)
            self.p += 8 * n
            return RObj("real", v, self.attrs(has_attr))
        if typ == 16:  # STRSXP
            n = self.length()
            v = [self.item().value for _ in range(n)]
            return RObj("str", v, self.attrs(has_attr))
        if typ in (19, 20):  # VECSXP, EXPRSXP
            n = self.length()
            v = [self.item() for _ in range(n)]
            return RObj("list", v, self.attrs(has_attr))
        if typ == 238:  # ALTREP
            info = self.item()
            state = self.item()
            attr = self.item()
            cls = info.value[0][1].value
            if cls == "compact_intseq":
                n, start, inc = state.value
                v = start + inc * np.arange(int(n))
                return RObj("int", v.astype(np.int64))
            if cls == "compact_realseq":
                n, start, inc = state.value
                return RObj("real", start + inc * np.arange(int(n)))
            if cls.startswith("wrap_"):
                return state.value[0]
            if cls == "deferred_string":
                arg = state.value[0][1]
                return RObj("str", [str(x) for x in arg.value])
            raise ValueError(f"ALTREP class {cls} not supported")
        raise ValueError(f"type {typ} not supported")


def read_rda(path):
    raw = open(path, "rb").read()
    for opener in (gzip.decompress, bz2.decompress, lzma.decompress):
        try:
            raw = opener(raw)
            break
        except Exception:
            pass
    assert raw[:5] in (b"RDX2\n", b"RDX3\n"), raw[:5]
    r = Reader(raw[5:])
    assert r.b[:2] == b"X\n"
    r.p = 2
    ver = r.int()
    r.int(); r.int()
    if ver == 3:
        n = r.int()
        r.p += n
    top = r.item()
    return {tag: obj for tag, obj in top.value}


def as_columns(df):
    """Data frame RObj -> dict of numpy arrays (factors -> label strings)."""
    names = df.attr["names"].value
    out = {}
    for nm, col in zip(names, df.value):
        if "levels" in col.attr:
            lev = col.attr["levels"].value
            out[nm] = np.array([None if c == NA_INT else lev[c - 1] for c in col.value], dtype=object)
        elif col.typ in ("int", "lgl"):
            v = col.value.astype(float)
            v[col.value == NA_INT] = np.nan
            out[nm] = v
        elif col.typ == "str":
            out[nm] = np.array(col.value, dtype=object)
        else:
            out[nm] = col.value
    return out
