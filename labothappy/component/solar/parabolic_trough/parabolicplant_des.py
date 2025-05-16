import __init__
import numpy as np
import component.solar.parabolic_trough.parabolictroughcollector as parabolictroughcollector
from geometries.solar.parabolictrough_geometry import PT_Collector_Geom

PT_geom = PT_Collector_Geom()
PT_geom.set_parameters("Soponova_MicroCSP")

