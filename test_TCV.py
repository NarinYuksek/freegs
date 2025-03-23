#!/usr/bin/env python

import freegs
import matplotlib.pyplot as plt

#########################################
# Create the machine, which specifies coil locations
# and equilibrium, specifying the domain to solve over
#23.03.2025, Narin Yüksek, do a test run of TCV equilibrium using equilibrium parameters from https://github.com/MrLou1976/PHY_329_Final_Project/blob/4c468f9ab0128b34ee7134530a844c8c0a918a23/TokamakModel.py#L106

tokamak = freegs.machine.TCV()

eq = freegs.Equilibrium(tokamak=tokamak,
                        Rmin=0.1, Rmax=2.0,    # Radial domain
                        Zmin=-2.0, Zmax=2.0,   # Height range
                        nx=129, ny=129)          # Number of grid points

#########################################
# Plasma profiles

profiles = freegs.jtor.ConstrainPaxisIp(eq,
                                        3e3, # Plasma pressure on axis [Pascals]
                                        -150000, # Plasma current [Amps]
                                        0.4) # vacuum f = R*Bt

#########################################
# Coil current constraints
#
# Specify locations of the X-points
# to use to constrain coil currents

xpoints = [(0.7, -1.1),   # (R,Z) locations of X-points
           (0.7, 1.1)]

isoflux = [(0.7, -1.1, 1.45, 0.0), (0.7, 1.1, 1.45, 0.0)] # (R1,Z1, R2,Z2) pair of locations

constrain = freegs.control.constrain(xpoints=xpoints, gamma=1e-12, isoflux=isoflux)

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

with open("diiid.geqdsk", "w") as f:
    geqdsk.write(eq, f)
