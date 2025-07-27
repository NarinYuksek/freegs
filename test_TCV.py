#!/usr/bin/env python

import freegs
import matplotlib.pyplot as plt
import numpy as np

#########################################
# Create the machine, which specifies coil locations
# and equilibrium, specifying the domain to solve over
#23.03.2025, Narin Yüksek, do a test run of TCV equilibrium using equilibrium parameters from https://github.com/MrLou1976/PHY_329_Final_Project/blob/4c468f9ab0128b34ee7134530a844c8c0a918a23/TokamakModel.py#L106

tokamak = freegs.machine.TCV()

eq = freegs.Equilibrium(tokamak=tokamak,
                        Rmin=0.1, Rmax=2.0,    # Radial domain
                        Zmin=-2.0, Zmax=2.0,   # Height range
                        nx=129, ny=129)        # Number of grid points

#########################################
# Plasma profiles

profiles = freegs.jtor.ConstrainPaxisIp(eq,
                                        3e3, # Plasma pressure on axis [Pascals] #originally 3e3
                                        -150000, # Plasma current [Amps] #originally 150000
                                        0.4) # vacuum f = R*Bt #originally 0.4

R_wall = [1.13600, 1.13600, 0.96790, 0.67070, 0.62400, 0.62400, 0.62400, 0.67240, 0.96790, 1.13600, 1.13600]
z_wall = [0.00000, 0.54940, 0.75000, 0.75000, 0.70330, 0.00000, -0.70330, -0.75000, -0.75000, -0.54940, 0.00000]

R_strike_lower = np.array([0.8])
R_strike_upper = R_strike_lower
z_strike_lower = np.array([-0.75])
z_strike_upper = -z_strike_lower

delta = -0.5
kappa = 1.5
R0 = 0.89
z0 = 0
r = 0.18

#Miller paramaterization of near last closed flux surface from https://doi.org/10.1063/1.872666

psis = 1
theta = np.linspace(0, 2*np.pi, 50)
Rs = R0 + psis*r*np.cos(theta + np.asin(delta)*np.sin(theta))
Zs = kappa*psis*r*np.sin(theta) + z0
psivals = []

for i in range (len(theta)):
    psivals.append((Rs[i], Zs[i], psis))

#From tcv_common.py
#INNER_LIMITER_R = 0.62400001
#OUTER_LIMITER_R = 1.14179182
#LIMITER_WIDTH = OUTER_LIMITER_R - INNER_LIMITER_R
#LIMITER_RADIUS = LIMITER_WIDTH / 2
#VESSEL_CENTER_R = INNER_LIMITER_R + LIMITER_RADIUS

#Coil current limits (min_current= -turns*ps_lim, and max_current= turns*ps_lim where N: number of turns, and ps_lim: power supply current limit)
current_lims= [(-7514, +7514), (-7514, +7514), (-7514, +7514), (-7514, +7514), (-7514, +7514), (-7514, +7514), (-7514, +7514), (-7514, +7514), (-7488, +7488), (-7488, +7488), (-7488, +7488), (-7488, +7488), (-7488, +7488), (-7488, +7488), (-7488, +7488), (-7488, +7488), (-7501, 7501), (-7501, 7501), (-7501, 7501), (-26026, +26026), (-26013, +26013), (-10764, +10764), (-7176, +7176)]

scale_fac= 10.0
for i in range(len(current_lims)):
    current_lims[i] = (current_lims[i][0]*scale_fac, current_lims[i][1]*scale_fac)

"""coils = [
        #             ("A1", Coil(0.422500,0.000000)),
        #             ("B1", Coil(0.445700,-0.936000)),
        #             ("B2", Coil(0.445700,0.936000)),
        #             ("C1", Coil(0.621500,-1.110000)),
        #             ("C2", Coil(0.621500,1.110000)),
        #             ("D1", Coil(1.176500,-1.170000)),
        #             ("D2", Coil(1.176500,1.170000)),
        ("E1", Coil(0.505000, -0.700000, turns=34)), #https://drive.google.com/file/d/1SOr3yDZZTjZ5mmAGGIBmAdx50DB7dxsl/view?usp=drive_link
        ("E2", Coil(0.505000, -0.500000, turns=34)),
        ("E3", Coil(0.505000, -0.300000, turns=34)),
        ("E4", Coil(0.505000, -0.100000, turns=34)),
        ("E5", Coil(0.505000, 0.100000, turns=34)),
        ("E6", Coil(0.505000, 0.300000, turns=34)),
        ("E7", Coil(0.505000, 0.500000, turns=34)),
        ("E8", Coil(0.505000, 0.700000, turns=34)),
        ("F1", Coil(1.309500, -0.770000, turns=36)),
        ("F2", Coil(1.309500, -0.610000, turns=36)),
        ("F3", Coil(1.309500, -0.310000, turns=36)),
        ("F4", Coil(1.309500, -0.150000, turns=36)),
        ("F5", Coil(1.309500, 0.150000, turns=36)),
        ("F6", Coil(1.309500, 0.310000, turns=36)),
        ("F7", Coil(1.309500, 0.610000, turns=36)),
        ("F8", Coil(1.309500, 0.770000, turns=36)),
        #             ("G1", Coil(1.098573,-0.651385)),
        #             ("G2", Coil(1.114000,-0.633000)),
        #             ("G3", Coil(1.129427,-0.614615)),
        #             ("G4", Coil(1.129427,0.614615)),
        #             ("G5", Coil(1.114000,0.633000)),
        #             ("G6", Coil(1.098573,0.651385)),
        ("T1", Coil(1.554000, -0.780000, turns=35)), #estimated turns, no info available
        ("T2", Coil(1.717000, -0.780000, turns=35)),
        ("T3", Coil(1.754000, -0.780000, turns=35)),
        ("OH1", Solenoid(0.43, -0.93, 0.93, 143)), #20.07.2025: Corrected number of turns (originally 100)
        (
            "OH2",
            Circuit(
                [
                    ("B1", Coil(0.445700,-0.936000, turns=29), 1.0),
                    ("B2", Coil(0.445700,0.936000, turns=29), 1.0)
                ]
            ),
        ),
        (
            "OH3",
            Circuit(
                [
                    ("C1", Coil(0.621500, -1.110000, turns=12), 1.0),
                    ("C2", Coil(0.621500, 1.110000, turns=12), 1.0)
                ]
            ),
        ),
        (
            "OH4",
            Circuit(
                [
                    ("D1", Coil(1.176500, -1.170000, turns=8), 1.0),
                    ("D2", Coil(1.176500, 1.170000, turns=8), 1.0)
                ]
            ),
        )
    ]"""

#########################################
# Coil current constraints
#
# Specify locations of the X-points
# to use to constrain coil currents

#xpoints = [(0.7, -1.1),   # (R,Z) locations of X-points
#          (0.7, 1.1)]  

#xpoints = [(0.7, -0.4),   # (R,Z) locations of X-points
#          (0.7, 0.4)]

#xpoints = [(1.0, -0.25),   # (R,Z) locations of X-points
#          (1.0, 0.25)]

xpoints = [(R0 + r*np.cos(np.pi/2 +  np.asin(delta)) + 0.00, kappa*r + 0.00 + z0), 
            (R0 + r*np.cos(np.pi/2 +  np.asin(delta)) + 0.00, -kappa*r - 0.00 + z0)]

#outer_mid = (1.09, 0.0)
#inner_mid = (0.75, 0.0)

outer_mid = [(R0 + r + 0.00), 0.0 + z0]
inner_mid = [(R0 - r - 0.00), 0.0 + z0]

#isoflux = [(0.7, -1.1, 1.45, 0.0), (0.7, 1.1, 1.45, 0.0)] # (R1,Z1, R2,Z2) pair of locations

#isoflux = [(0.7, -0.4, 1.0, 0.0), (0.7, 0.4, 1.0, 0.0)] # (R1,Z1, R2,Z2) pair of locations

isoflux = [(xpoints[0][0], xpoints[0][1], outer_mid[0], outer_mid[1]), (xpoints[1][0], xpoints[1][1], outer_mid[0], outer_mid[1]), (xpoints[0][0], xpoints[0][1], inner_mid[0], inner_mid[1]),
            (xpoints[0][0], xpoints[0][1], R_strike_lower[0], z_strike_lower[0]), (xpoints[1][0], xpoints[1][1], R_strike_upper[0], z_strike_upper[0])] #(0.82, -0.62, outer_mid[0], outer_mid[1])] # (R1,Z1, R2,Z2) pair of locations

#for i in range (1, len(theta)):
#    isoflux.append((psivals[0][0], psivals[0][1], psivals[i][0], psivals[i][1]))

#constrain = freegs.control.constrain(xpoints=xpoints, gamma=1e-12, isoflux=isoflux, psivals=psivals)

constrain = freegs.control.constrain(xpoints=xpoints, gamma=1e-12, isoflux=isoflux)

#constrain = freegs.control.constrain(xpoints=xpoints, gamma=1e-12, isoflux=isoflux, current_lims=current_lims)

#constrain = freegs.control.constrain(xpoints=xpoints, gamma=1e-12, current_lims=current_lims)

constrain(eq)

#########################################
# Nonlinear solve

freegs.solve(eq,          # The equilibrium to adjust
             profiles,    # The plasma profiles
             constrain,   # Plasma control constraints
             show=True)   # Shows results at each nonlinear iteration

# eq now contains the solution

print("Done!")

print("Plasma current: %e Amps" % (eq.plasmaCurrent()))
print("Pressure on axis: %e Pascals" % (eq.pressure(0.0)))
print("Plasma poloidal beta: %e" % (eq.poloidalBeta()))
print("Plasma volume: %e m^3" % (eq.plasmaVolume()))

eq.tokamak.printCurrents()

# plot equilibrium
axis = eq.plot(show=False)
axis.plot(R_wall, z_wall) #"c", (0.85, 0.325, 0.098)
axis.plot(Rs, Zs)
tokamak.plot(axis=axis, show=False)
constrain.plot(axis=axis, show=True)

# Safety factor
plt.plot(*eq.q())
plt.xlabel(r"Normalised $\psi$")
plt.ylabel("Safety factor")
plt.grid()
plt.show()

##############################################
# Save to geqdsk file

from freegs import geqdsk

with open("tcv.geqdsk", "w") as f:
    geqdsk.write(eq, f)
