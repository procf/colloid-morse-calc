# colloid-morse-calc
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/procf/colloid-morse-calc/HEAD)

## What is it?
**colloid-morse-calc** is a Python script for calculating the Morse Potential from experimental values. <br>
Have you ever wanted to replicate a specific colloidal depletion gel in simulation? This can help! <br>
See more info about assumptions and limitations below.

## Documentation
This code uses experimental details about a depletion-based spherical-colloid system to calculate the effective Morse Potential and related MD simulation parameters.

**A standard simulation with Morse Potential assumes that your particles are stabilized to behave like hard spheres.** They are not sticky or patchy, they do not interact with each other outside of depletion, their surface charge (including charge from any surface coating or polymer brush) is balanced by added salt, etc.

You can run the code in a Jupyter Notebook on Binder, but no data is saved between Binder sessions. <br>
Copy the output into another file, or fork and clone this repository to save your parameters and make changes to the code.

- [what you need](/README.md#what-you-need) (inputs)<br>
- [what you get](/README.md#what-you-get) (outputs)<br>
- *Jump to Morse Calculator*<br>
- *What if I don't know $c^*$ for this experimental system?*<br>
- *Extensions*<br>
- *More about the Morse Potential*<br>
- *Mathematical Form*<br>
- *More about simulations*<br>

### What you need
You need to provide the following information and **experiment** and **simulation**:

| Experimental value    | units   |
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

| Simulation value    | units   |
| :-------- | :------- |
| simulation timestep   | $dt$ |
| simulation recording interval | $trigger$ = $period$ |
| average colloid radii in simulation units | $R_{C,sim}$ |
| simulation solvent viscosity parameter | $eta_{0,sim}$ |
| simulation solvent number density | $\rho_{S,sim}$ |
| simulation temperature | $kT$ |


## What you get 
The calculator ouputs two sets of information: **real values** and **simulation values**

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
| Gravitational force on colloids | $F_G$ *[J = kg m<sup>2</sup>]* |

| Simulation value    | units   |
| :-------- | :------- |
| Morse potential attraction strength (rounded Depletion minimum)  | $D_0 = U_0/kT$ |
| Morse potential attraction range parameter | $\alpha = \kappa$ |
| real-time per simulation timestep | $dt$ [s] |
| real-time per recorded simulation frame | $t_{frame}$ [s] |
| number of frames needed to match experiment time | $nframes$ |
| simulation unit length | $[length]$ [m] |
| simulation unit mass | $[mass]$ [g] |
| simulation unit energy | $[energy]$ [J] |
| simulation unit time | $[time]$ [s] |
| standard gravitational acceleration in simulation units | $g_{sim}$ |
| Gravitational Peclet number in simulation units | $Pe_{G,sim}$ |
| force of gravity per unit length, in sim units | $F_{G,sim}$ |
