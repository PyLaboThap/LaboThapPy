# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 14:09:18 2026

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

def plot_blank_Ts(fluid, xlim=None, ylim=None):
    """
    Create a T-s diagram with saturation curve for the given fluid.
    Returns the matplotlib figure and axes without showing the plot.
    
    Parameters
    ----------
    fluid : str
        Name of the fluid as recognized by CoolProp (e.g., "Water", "R134a").
    xlim : tuple of float, optional
        Limits for the x-axis (entropy), e.g., (s_min, s_max)
    ylim : tuple of float, optional
        Limits for the y-axis (temperature), e.g., (T_min, T_max)
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    ax : matplotlib.axes.Axes
        The axes object.
    """
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Set labels and title
    ax.set_xlabel("Entropy s [J/kg/K]")
    ax.set_ylabel("Temperature T [K]")
    ax.set_title(f"T-s Diagram with Saturation Curve for {fluid}")
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Get triple point and critical temperatures
    T_triple = CP.PropsSI('T_triple', fluid)
    T_crit   = CP.PropsSI('T_critical', fluid)
    
    # Temperature array from triple point up to just below critical
    T_vals = np.linspace(T_triple, T_crit * 0.99999, 500)
    
    # Saturated liquid (Q=0) and vapor (Q=1) entropy
    s_liq = [CP.PropsSI('S', 'T', T, 'Q', 0, fluid) for T in T_vals]
    s_vap = [CP.PropsSI('S', 'T', T, 'Q', 1, fluid) for T in T_vals]
    
    # Plot saturation curves
    ax.plot(s_liq, T_vals, color="black")
    ax.plot(s_vap, T_vals, color="black")
    
    # Apply axis limits if provided
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    
    # Optional: legend
    ax.legend()
    
    return fig, ax


fig, ax = plot_blank_Ts('Propane', xlim = [-2000, 6000], ylim = [50, 500])

plt.show()
