def zero_brent(a, b, machep, t, f, tol_f):
    """
    Seeks the root of a function f(x) in an interval [a, b] using Brent's method.

    The interval [a, b] must be a change-of-sign interval for f:
    f(a) and f(b) must have opposite signs. Assuming f is continuous,
    this guarantees at least one root c in (a, b) with f(c) = 0.

    Parameters
    ----------
    a, b : float
        Endpoints of the change-of-sign interval.
    machep : float
        Estimate of relative machine precision.
    t : float
        Positive error tolerance.
    f : callable
        Function whose zero is being sought, f(x) -> float.
    tol_f : float
        Error tolerance targeted on f (addon by Remi Dickes, 26/11/2018).

    Returns
    -------
    value : float
        Estimated zero of f.
    f_out : float
        Value of f at the estimated zero.
    """

    # Make local copies of a and b
    sa = a
    sb = b
    fa = f(sa)

    if fa == 0 or abs(fa) < tol_f:
        return sa, fa

    fb = f(sb)
    if fb == 0 or abs(fb) < tol_f:
        return sb, fb

    c = sa
    fc = fa
    e = sb - sa
    d = e

    while True:

        if abs(fc) < abs(fb):
            sa = sb
            sb = c
            c = sa
            fa = fb
            fb = fc
            fc = fa

        tol = 2.0 * machep * abs(sb) + t
        m = 0.5 * (c - sb)

        if abs(m) <= tol or fb == 0.0 or abs(fb) < tol_f:
            break

        if abs(e) < tol or abs(fa) <= abs(fb):
            e = m
            d = e
        else:
            s = fb / fa

            if sa == c:
                p = 2.0 * m * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)

            if p > 0.0:
                q = -q
            else:
                p = -p

            s = e
            e = d

            if 2.0 * p < 3.0 * m * q - abs(tol * q) and p < abs(0.5 * s * q):
                d = p / q
            else:
                e = m
                d = e

        sa = sb
        fa = fb

        if tol < abs(d):
            sb = sb + d
        elif m > 0.0:
            sb = sb + tol
        else:
            sb = sb - tol

        fb = f(sb)

        if (fb > 0.0 and fc > 0.0) or (fb <= 0.0 and fc <= 0.0):
            c = sa
            fc = fa
            e = sb - sa
            d = e

    value = sb
    f_out = fb
    return value, f_out
