#!/usr/bin/python3

# --- loading libraries 

from CoolProp.CoolProp import PropsSI
from numpy import *
import matplotlib.pyplot as plt
import numpy as np

# --- iteration criteria 

tol = 1e-4
maxIter = 20

# --- constants

radToDeg = 180/arccos(-1)

# --- cycle parameters 

fluid = 'Cyclopentane'
mdot  = 34.51        # mass flow rate [kg/s]
p01   = 767800       # inlet total pressure [Pa]
T01   = 273.15 + 131 # inlet total temperature [K]
pn    = 82000          # outlet pressure [Pa]

P = 2506000 # work to be performed 

# --- design parameters 

workCoeff = 1               # ideal point in the Schmidt diagram 
flowCoeff = 0.6             # " 
reaction  = 1 - workCoeff/2 # axial leaving stage
Mmax = 0.8                  # maximum Mach number allowed between the blade rows

"""
Blade row efficiency adjusted to get the right pressure ratio
Higher than researched turbine total-to-static efficiency because of the exit losses
"""

bladeRowEfficiency = 0.818  # blade row / static-to-static efficiency (h1-h2)/(h1-h2s)
Zweifel = 0.8               # choose moderate loading
ReMin = 5e5                 # turbulent flow
minAspectRatio = 1          # minimum ratio blade height / blade chord

# ==============================================================================
# --- turbomachinery computation routines
# ==============================================================================

# --- compute velocity triangle shape in function of stage coefficients --------
# --- velocity components related to the blade speed u -------------------------
def computeVelTriangle(psi,phi,R):
# ------------------------------------------------------------------------------

    vu2OverU = (2*(1-R)+psi)/2
    vu3OverU = (2*(1-R)-psi)/2
    vmOverU  = phi
    return (vu2OverU,vu3OverU,vmOverU)

# --- compute a single blade passage given (relative) velocities up and down ---
def computeBladeRow(h1,s1,cm,cu1,cu2,eff):
# ------------------------------------------------------------------------------
    
    h01 = h1  + (cu1*cu1 + cm*cm)/2
    h2  = h01 - (cu2*cu2 + cm*cm)/2
    h2s = h1  - (h1-h2)/eff
    p2 = PropsSI('P','H',h2s,'S',s1,fluid)
    s2 = PropsSI('S','H',h2 ,'P',p2,fluid)
    
    return (h2,s2)

# --- compute static conditions in stage given velocity triangles --------------
def computeStage(h1,s1,vm,vu1,vu2,vu3,u,eff):
# ------------------------------------------------------------------------------

    wu2 = vu2 - u
    wu3 = vu3 - u
    
    (h2,s2) = computeBladeRow(h1,s1,vm,vu1,vu2,eff)     # --- stator row 
    (h3,s3) = computeBladeRow(h2,s2,vm,wu2,wu3,eff) # --- rotor row 

    return (h2,s2,h3,s3)

# --- compute n stages ---------------------------------------------------------
def computeRepeatingStages(h1,s1,vm,vu1,vu2,vu3,u,eff,n):
# ------------------------------------------------------------------------------
    
    stages = []

    hn = h1
    sn = s1
    
    for i in range(int(n)):

        sc = computeStage(hn,sn,vm,vu1,vu2,vu3,u,eff)
        stages.append(sc)
        
        print("Stage %i : (h1,s1) = (%g,%g) -> (h2,s2) = (%g,%g) -> (h3,s3) = (%g,%g)"
              %(i,hn,sn,sc[0],sc[1],sc[2],sc[3]))
        
        hn = sc[2]
        sn = sc[3]
        
    return stages

# ==============================================================================
# --- Compute the work required and the estimated total-to-static efficiency
# ==============================================================================

# --- compute the inlet total conditions

h01 = PropsSI("H","P",p01,"T",T01,fluid)
s01 = PropsSI("S","P",p01,"T",T01,fluid)
c01 = PropsSI("A","H",h01,"S",s01,fluid)
d01 = PropsSI('D','P',p01,'T',T01,fluid)

print("Inlet total enthalpy h01 = %g kJ/kg, entropy s01 = %g kJ/kg.K, d01 = %g kg/m3, a01 = %g m/s"
      %(h01/1000,s01/1000,d01,c01))

# --- isentropic total to static expansion

ss = s01
hs = PropsSI("H","S",ss,"P",pn,fluid)
ps = PropsSI('P','H',hs,'P',pn,fluid)
Ts = PropsSI('T','H',hs,'P',pn,fluid)
cs = PropsSI('A','H',hs,'P',pn,fluid)
ds = PropsSI('D','H',hs,'P',pn,fluid)

print("Isentropic expansion ps = %g Pa, Ts = %g K, ds = %g kg/m3, hs = %g kJ/kg, ss = %g kJ/kg.K, as = %g m/s" %(ps,Ts,ds,hs/1e3,ss/1e3,cs))

dh0s = h01 - hs

print("Isentropic work %g [kJ/kg] and power %g [MW], p01/pn = %g"%(dh0s/1e3,mdot*dh0s/1e6,p01/pn))

# --- total-to-total conditions for the real expansion

dh0 = P/mdot
eff = dh0/dh0s
h0n = h01 - dh0
print("Isentropic total-to-static efficiency targeted %g and work dh0 = %g kJ/kg"%(eff,dh0/1e3))

# ==============================================================================
# --- Compute velocity triangle shape for given stage coefficients
# ==============================================================================

(vu2,vu3,vm) = computeVelTriangle(workCoeff,flowCoeff,reaction)

wu2 = vu2 - 1
wu3 = vu3 - 1

v2 = sqrt(vu2*vu2+vm*vm)
w2 = sqrt(wu2*wu2+vm*vm)
v3 = sqrt(vu3*vu3+vm*vm)
w3 = sqrt(wu3*wu3+vm*vm)

a2 = arctan(vu2/vm)
b2 = arctan(wu2/vm)
a1 = a3 = arctan(vu3/vm)
b1 = b3 = arctan(wu3/vm)

solidityStator = 2*cos(a2)/cos(a1)*sin(abs(a2-a1))/Zweifel
solidityRotor  = 2*cos(b3)/cos(b2)*sin(abs(b3-b2))/Zweifel

print("Velocity triangle shape:")
print("vu2/u = %g, wu2 = %g, vu3/u = %g, wu3 = %g, vm/u = %g, v2/u = %g, w3/u = %g"
      %(vu2,wu2,vu3,wu3,vm,v2,w3))
print("a2 = %g°, b2 = %g°, a3 = %g°, b3 = %g°"
      %(a2*radToDeg,b2*radToDeg,a3*radToDeg,b3*radToDeg))
print("Solidity stator = %g and rotor %g"%(solidityStator,solidityRotor))

# ==============================================================================
# determine rotation speed by limiting maximum Mach number it is found
# either at stator or rotor outlet apparently, the speed of sound is
# pretty much constant so for simplicity we take the one corresponding
# to the isentropic expansion
# ==============================================================================

vMax = cs * Mmax
u = vMax / max([v2,w3])
dh0Stage = workCoeff * u * u
nStages = round(dh0/dh0Stage)
dh0Stage = dh0/nStages
u = sqrt(dh0Stage/workCoeff)

print("Rotation speed u = %g m/s, dh0 = %g kJ/kg, number of stages %g"
      %(u,dh0Stage/1000,dh0/dh0Stage))

vm  *= u
vu2 *= u
vu3 *= u
wu2 *= u
wu3 *= u
vu1 = vu3

exit_loss = (vm*vm+vu3*vu3)/2

print("Velocities: "
      "vm = %g m/s, vu2 = %g m/s, wu2 = %g m/s, vu3 = %g m/s, wu3 = %g m/s"
      %(vm,vu2,wu2,vu3,wu3))

print("Exit loss = %g kJ/kg"%(exit_loss/1e3))


# ==============================================================================
# --- Compute the exit conditions from the estimated conditions and power
# ==============================================================================

h1 = h01 - vm*vm/2
s1 = s01
p1 = PropsSI("P","H",h1,"S",s1,fluid)
T1 = PropsSI("T","H",h1,"S",s1,fluid)
c1 = PropsSI("A","H",h1,"S",s1,fluid)
d1 = PropsSI('D','P',p1,'T',T1,fluid)

print("Isentropic expansion "
      "p1 = %g P, T1 = %g K, d1 = %g kg/m3, h1 = %g kJ/kg, s1 = %g kJ/kg.K, a1 = %g m/s"
      %(p1,T1,d1,h1/1e3,s1/1e3,c1))

# --- real expansion total conditions 

hn  = h0n - vm*vm/2
Tn = PropsSI('T','H',hn,'P',pn,fluid)
sn = PropsSI('S','H',hn,'P',pn,fluid)
cn = PropsSI('A','H',hn,'P',pn,fluid)
dn = PropsSI('D','H',hn,'P',pn,fluid)

print("Real expansion static conditions "
      "pn = %g Pa, Tn = %g K, dn = %g kg/m3, hn = %g kJ/kg, sn = %g kJ/kg.K, an = %g m/s"
      %(pn,Tn,dn,hn/1e3,sn/1e3,cn))

s0n = sn
T0n = PropsSI('T','H',h0n,'S',s0n,fluid)
p0n = PropsSI('P','H',h0n,'S',s0n,fluid)
c0n = PropsSI('A','H',h0n,'S',s0n,fluid)
d0n = PropsSI('D','H',h0n,'S',s0n,fluid)

print("Real expansion total conditions "
      "p0n = %g Pa, T0n = %g K, d0n = %g kg/m3, h0n = %g kJ/kg, s0n = %g kJ/kg.K, a0n = %g m/s"
      %(p0n,T0n,d0n,h0n/1e3,s0n/1e3,c0n))


# ==============================================================================
# --- Compute conditions through stages given work and efficiency
# ==============================================================================

h_s = computeRepeatingStages(h1,s1,vm,vu1,vu2,vu3,u,bladeRowEfficiency,nStages)

# 1) blade row efficiency to get to the pressure ratio for fixed work
# Calcul de Pr_ts : fitter avec eta_aubage

p_out_static = PropsSI('P','H',h_s[-1][-2], 'S', h_s[-1][-1], fluid)
res = p_out_static - pn # Iterate on blade efficiency

# 2) Determine densities and throughflow area

nStages = int(nStages)

rho_out = np.zeros([nStages,2])
A_flow = np.zeros([nStages,2])

mu_out = np.zeros([nStages,2])
cord_min = np.zeros([nStages,2])
cord_solid = np.zeros([nStages,2])

pitch = np.zeros([nStages,2])
n_blade = np.zeros([nStages,2])

h_blade = np.zeros([nStages,2])
AR = np.zeros([nStages,2])

r_moy = 0.2 # Iterate on this value so that AR > 1

AR_imposed = np.array([[1, 1.21],
                       [1.43, 1.64],
                       [1.86, 2.07],
                       [2.29, 2.5]])

for i in range(int(nStages)):
      rho_out[i][0] = PropsSI('D','H',h_s[i][0], 'S', h_s[i][1], fluid)
      rho_out[i][1] = PropsSI('D','H',h_s[i][2], 'S', h_s[i][3], fluid)

      A_flow[i][0] = mdot/(rho_out[i][0]*vm)
      A_flow[i][1] = mdot/(rho_out[i][1]*vm)

# 3) determine minimum chord to satisfy minimum Reynolds
#    by using velocity, density and viscosity at the blade outlet

      mu_out[i][0] = PropsSI('V','H',h_s[i][0], 'S', h_s[i][1], fluid)
      mu_out[i][1] = PropsSI('V','H',h_s[i][2], 'S', h_s[i][3], fluid)

      h_blade[i][0] = A_flow[i][0]/(4*np.pi*r_moy)
      h_blade[i][1] = A_flow[i][1]/(4*np.pi*r_moy)

      cord_min[i][0] = (ReMin*mu_out[i][0])/(rho_out[i][0]*vm)
      cord_min[i][1] = (ReMin*mu_out[i][1])/(rho_out[i][1]*vm)

      AR_max = min(h_blade/cord_min) 

AR_imposed = np.linspace(1,AR_max,nStages*2).reshape(4,2)
cord_calc = h_blade/AR_imposed

for i in range(int(nStages)):

      pitch[i][0] = solidityStator*cord_calc[i][0]
      pitch[i][1] = solidityRotor*cord_calc[i][1]

      n_blade[i][0] = round(2*np.pi*r_moy/pitch[i][0])
      n_blade[i][1] = round(2*np.pi*r_moy/pitch[i][1])

      # A_flow = pi*((r_m + h)**2 - (r_m - h)**2)
      #        = pi*((r_m**2 - r_m**2) + (2*r_m*h - (-2*r_m*h)) + (h**2 - h**2))
      #        = pi*4*r_m*h
      # h = A_flow/(4*pi*r_m)

      AR[i][0] = h_blade[i][0]/cord_calc[i][0]
      AR[i][1] = h_blade[i][1]/cord_calc[i][1]



# 5) determine rotation speed
omega = u/r_moy
omega_RPM = omega*60/(2*np.pi)

# 6) Compute geometry

r_hub = []
r_ext = []
r_ratio2 = []
r_hub_tip = []

for i in range(nStages):    
      r_ext.append(r_moy + h_blade[i][0]/2)
      r_hub.append(r_moy - h_blade[i][0]/2)
      r_ratio2.append((r_ext[-1]/r_hub[-1])**2)
      r_hub_tip.append(r_hub[-1]/r_ext[-1])

      r_ext.append(r_moy + h_blade[i][1]/2)
      r_hub.append(r_moy - h_blade[i][1]/2)
      r_ratio2.append((r_ext[-1]/r_hub[-1])**2)
      r_hub_tip.append(r_hub[-1]/r_ext[-1])

# 7) Compute Thermodynamic properties

p_out = [p1]
s_out = [s1]
h_out = [h1]

for i in range(int(nStages)):
      p_out.append(PropsSI('P','H',h_s[i][0], 'S', h_s[i][1], fluid))
      p_out.append(PropsSI('P','H',h_s[i][2], 'S', h_s[i][3], fluid))

      s_out.append(h_s[i][1])
      s_out.append(h_s[i][3])

      h_out.append(h_s[i][0])
      h_out.append(h_s[i][2])

# 8) Plot Results

r_moy_line = np.ones(len(r_ext))*r_moy

x = np.linspace(0,len(r_ext)-1, len(r_ext))

labels = []
i = 1

while len(labels) < len(x):
      labels.append("S" + str(i))
      labels.append("R" + str(i))
      i += 1

fontsize = 16
ticksize = 12


plt.figure()
plt.plot(r_ext)
plt.plot(r_hub)
plt.plot(r_moy_line)
plt.plot(h_blade.flatten())
#  plt.axis([-0.5, len(r_ext)-0.5, 0, 0.3])
plt.legend(["$r_{ext}$", "$r_{hub}$", "$r_{moy}$", "$h_{blade}$"])
plt.xticks(ticks=x, labels=labels, size=ticksize)
plt.grid(axis = 'x')
plt.ylabel("Length or radius [m]", fontsize= fontsize)
plt.show()

plt.figure()
plt.plot(AR.flatten())
# plt.axis([-0.5, len(r_ext)-0.5, 0, 8])
plt.xticks(ticks=x, labels=labels, size=ticksize)
plt.grid()
plt.ylabel("Aspect Ratio [-]", fontsize= fontsize)
plt.show()

n_blade_plot = n_blade.flatten()

for i in range(len(n_blade_plot)*2):
      if mod(i,4) == 1 or mod(i,4) == 2: # Stator
            n_blade_plot = np.insert(n_blade_plot,i,None)      

n_blade_plot = n_blade_plot.reshape(int(len(n_blade_plot)/2),2)

plt.figure()
plt.plot(n_blade_plot[:, 0], 'o', label="Stator Blades")  # Plot first column with label
plt.plot(n_blade_plot[:, 1], 'o', label="Rotor Blades")  # Plot second column with label
# plt.axis([-0.5, len(r_ext)-0.5, 0, 8])
plt.xticks(ticks=x, labels=labels, size=ticksize)
plt.legend()
plt.grid()
plt.ylabel("Blade number [-]", fontsize= fontsize)
plt.show()

plt.figure()
plt.plot(r_ratio2)
# plt.axis([-0.5, len(r_ext)-0.5, 0, 30])
plt.xticks(ticks=x, labels=labels, size=ticksize)
plt.grid()
plt.ylabel("$\\left[ r_{ext}/r_{hub} \\right]^2$ [-]", fontsize= fontsize)
plt.show()

plt.figure()
plt.plot(r_hub_tip)
# plt.axis([-0.5, len(r_ext)-0.5, 0, 30])
plt.xticks(ticks=x, labels=labels, size=ticksize)
plt.grid()
plt.ylabel("$\\left[ r_{hub}/r_{tip} \\right]$ [-]", fontsize= fontsize)
plt.show()

# Thermo Prop
x2 = np.linspace(0,len(r_ext), len(r_ext)+1)
labels2 = ['0'] + labels

plt.figure()
plt.plot(np.array(p_out)*1e-3)
# plt.axis([-0.5, len(r_ext)+0.5, 0, 1*1e3])
plt.xticks(ticks=x2, labels=labels2, size=ticksize)
plt.grid()
plt.ylabel("$Pressure$ [kPa]", fontsize= fontsize)
plt.show()

plt.figure()
plt.plot(s_out, h_out)
plt.plot([s1, s1], [h1, hs])

# Define entropy range (in J/(kg·K))
entropy_range = np.linspace(s_out[0], s_out[-1], 100)  # Adjust range for your fluid

for P in p_out:
    enthalpy = [PropsSI('H', 'S', s, 'P', P, fluid) for s in entropy_range]  # Enthalpy in kJ/kg
    entropy = entropy_range  # Entropy in kJ/(kg·K)
    plt.plot(entropy, enthalpy, color = 'grey', alpha=0.3, label=f'P = {P/1e5} bar')

plt.ylabel("$Enthalpy$ [J/kg]", fontsize= fontsize)
plt.xlabel("$Entropy$ [J/(kg x K)]", fontsize= fontsize)
plt.legend(["real", "isentropic"])
plt.show()


print(f"Turbine mean diameter: {r_moy} [m]")
print(f"Turbine rotation speed: {omega_RPM} [RPM]")
print(f"Turbine number of stage : {nStages} [-]")
print(f"Turbine static-to-static blade efficiency : {bladeRowEfficiency} [-]")
