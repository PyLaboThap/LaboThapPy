def stepsmoother(x0, x1, x):
    """
    Smooth cubic Hermite interpolation between two regimes.
    
    Returns 0 for x <= x0, 1 for x >= x1, and a smooth cubic blend in between.
    The function is continuous with zero derivatives at boundaries.

    This is used to smoothly transition between laminar, transition, and turbulent
    flow regimes without discontinuities.

    Parameters
    ----------
    x0 : float
        Lower boundary of interpolation interval [-]
    x1 : float
        Upper boundary of interpolation interval [-]
    x : float
        Input value to interpolate [-]

    Returns
    -------
    float
        Interpolated weight: 0 <= results <= 1    
    """
    if x <= x0:
        return 0.0
    elif x >= x1:
        return 1.0
    else:
        t = (x - x0) / (x1 - x0)
        return t**2 * (3 - 2*t)