# colloid-morse-calc
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/procf/colloid-morse-calc/HEAD)

## What is it?
**colloid-morse-calc** is an interactive Python script for calculating the Morse Potential from experimental values. <br>
Have you ever wanted to replicate a specific colloidal depletion gel in simulation? This can help! <br>
See more info about assumptions and limitations below.

## Documentation
This code uses information about a depletion-based spherical-colloid experiment to calculate a set of parameters for the Morse Potential and related MD simulation values.

**A standard simulation with Morse Potential assumes that your particles are stabilized to behave like hard spheres.** They are not sticky or patchy, they do not interact with each other outside of depletion, their surface charge (including charge from any surface coating or polymer brush) is balanced by added salt, etc.

*What you need*<br>
*What you get*<br>
*Jump to Calculation*<br>
*What if I don't know $c^*$ for this experimental system?*<br>
*More about the Morse Potential*<br>
*Mathematical Form*<br>
*More about simulations*<br>

#### You need to provide the following information:

| value    | units   |
| :-------- | :------- |
| depletant molar mass   | $M_d$ *[g/Mol]* |
| depletant concentration | $c_d$ *[mg/mL = g/L]* |
| depletant overlap concentration     | $c^*$ *[mg/mL = g/L]* |
| solvent density  | $\rho_s$ *[g/mL = kg/L]* |
| solvent dynamic viscosity  | $\eta$ *[Pa s = (kg/(m s<sup>2</sup>)) s]* |
| colloid particle volume fraction  | $\phi$ |
| average colloid particle diameter(s)  | $d_c$ *[meter]* |
| colloid particle density | $\rho_c$ *[g/mL = kg/L]* |
| time (between sample mixing and data capture)  | $t$ *[sec]* |

#### What you get: Real values, and Simulation Values

| Real value    | units   |
| :-------- | :------- |
| depletant radius of gyration   | $r_g$ *[nm]* |
| depletant attraction range | $d_g = 2r_g$ *[nm]* |
| depletant-to-colloid size ratio  | $\Delta$ |
| Depletion potential minimum  | $D_0 = U_0/kT$ |
| bare colloid diffusion time  | $\tau_D$ *[s]* |
| data collection time in diffusion times | $t/\tau_D$ |
| approximate colloid mass(es)  | $m_C$ *[nanograms]* |
| Gravitational Peclet number [advection/settling] | $Pe_G$ |
| Gravitational force per colloid  | $F_G$ *[J = kg m<sup>2</sup>]* |

| Simulation value    | units   |
| :-------- | :------- |
| simulation unit length | $[length]$ [m] |
| simulation unit mass | $[mass]$ [g] |
| simulation unit energy | $[energy]$ [J] |
| simulation unit time | $[time]$ [s] |
| Morse potential attraction strength (rounded Depletion minimum)  | $D_0 = U_0/kT$ |
| Morse potential attraction range parameter | $\alpha = \kappa$ |
| real-time per simulation timestep | $dt$ [s] |
| real-time per recorded simulation frame | $t_{frame}$ [s] |
| number of frames needed to match experiment time | $nframes$ |
| standard gravitational acceleration in simulation units | $g_{sim}$ |
| Gravitational Peclet number in simulation units | $Pe_{G,sim}$ |
| force of gravity per unit length, in sim units | $F_{G,sim}$ |
