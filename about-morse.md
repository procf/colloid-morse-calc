### More about the Morse Potential
- The Morse Potential was introduced by Philip Morse in 1929 as an exact solution for harmonic interactions in a diatomic molecule</br>
  [DOI: 10.1103/PhysRev.34.57](https://doi.org/10.1103/PhysRev.34.57)
- In 2010, Taffs et.al. showed that the Morse Potential can describe an Asakura-Oosawa-Vrij (AO) type fluid with more accuracy than the standard one-component AO description used by the Depletion Potential<br>
  [DOI: 10.1088/0953-8984/22/10/104119](https://doi.org/10.1088/0953-8984/22/10/104119)
- In 2014, Zia et.al used the Morse potential to study the micro-mechanical coarsening and rheology of colloidal gels, establishing it as a model of choice for colloidal gel rheology<br>
  [DOI 10.1122/1.4892115](https://doi.org/10.1122/1.4892115)

**Assumptions**
- The Morse potential assumes that particles are hard spheres
- The depletant-to-colloid size ratio is assumed to be $\Delta<0.154$ 
  - The standard (one-component) AO Depletion model can only ignore higher order (3-body+) interactions when the ratio of depltant-to-colloid size is $\Delta=r_g/R_C<0.154$ ; for $\Delta=0.15$ Morse is an exact match for Depletion. Above this value there may be non-negligible many-body effects
- The attraction range is set by $\alpha=3/r_g$
  - The width of the Morse potential well is set by an $\alpha$ parameter; Zia et.al. showed that when $\Delta=0.1$, $\alpha R_C=30$ ; therefore, you can approximate $\alpha$ for any system using $\alpha=3/\Delta R_C=3/r_g$

<br>

### Mathematical Form: Morse vs. Depletion

| <div style="width:350px">Morse Potential</div> | <div style="width:350px">Depletion Potential (one-component AO)</div> |
|:-----------------:|:--------------------------:|
| $U_{morse}(r) = U_0 \left( e^{-2\alpha(r-r_0)} - 2e^{-\alpha(r-r_0)} \right) $ |  $U_{AO}(r)=-\Pi V_{ov}(r)$ |
| $U_0=U_{morse}(r_0)=D_0 k_B T$ <br> $r_0 = R_{Ci}+R_{Cj}$ <br> $\alpha = 3/r_g = 3/\Delta \langle R_C \rangle$ | $\Pi=n_bk_BT$ <br> $n_b=\left( \frac{1}{1-\phi_C}\right) \left( \frac{c_d N_A}{M} \right)$ </br> $V_{ov}(r)=\frac{\pi}{6} \left( 2r_g - h \right)^2 \left( 3\langle R_C \rangle +2r_g \frac{h}{2} \right)$ <br> $h = r-r_0$ <br> $\langle R_C \rangle = 2 \left( \frac{R_{Ci}R_{Cj}}{R_{Ci}+R_{Cj}} \right)$ <br> $r_g = \left( \frac{3M}{4\pi N_A c^*} \right)^{1/3}$|


| <div style="width:390px">Morse Parameters</div> | <div style="width:40px">symbol</div> |
|------------------|--------|
| attraction strength | $D_0$ | 
| attraction range parameter | $\alpha$ | 
| center-center particle separation distance | $r$ | 
| particle radii | $R_{Ci}$ ; $R_{Cj}$ |
| center-center distance at particle contact | $r_0 = R_{Ci} + R_{Cj}$ |
| energy scale (Boltzmann constant * Temperature) | $k_B T$ |

| <div style="width:400px">Depletion Parameters</div> | <div style="width:40px">symbol</div> |
|----------------------|--------|
| osmotic pressure (Van 't Hoff equation) | $\Pi$ |
| depletion zone overlap volume | $V_{ov}$ |
| depletant concentration | $c_d$|
| depletant overlap concentration | $c^*$|
| depletant size (radius of gyration) | $r_g$ |
| depletant number density | $n_b$ |
| colloidal volume fraction | $\phi_C$ |
| depletant molar mass: number average ($M_n$) or weight average ($M_w$) | $M$ |
| Avogadro's number | $N_A$ |
| depletant size (radius of gyration) | $r_g$ |
| center-center particle separation distance | $r$ |
| particle radii | $R_{Ci}$ ; $R_{Cj}$ |
| center-center distance at particle contact | $r_0 = R_{Ci} + R_{Cj}$ |
| surface-surface separation distance | $h=r-r_0$ |
| energy scale (Boltzmann constant * Temperature) | $k_B T$ |

<br>

### Simulating Morse
In simulation all parameters are non-dimensionalized, this calculation uses the following conversions:
  | parameter | units |
  |:-----------:|:-------:|
  | length | smallest colloidal particle radius $R_{C,min}$|
  | mass | colloid mass $m_C$ *(see note 1) |
  | energy | $k_B T = 0.1$ |
  | time | bare colloid diffusion time *(see notes 2,3)|

*Note 1:* for density matching with solvent in DPD simulations, mass units are also scaled by colloid volume times a system-wide "solvent" number-density value: $V_C \cdot \rho_{sim}$ <br>

*Note 2:* Bare colloid diffusion time can be measured in experiment or approximated from the Stokes–Einstein–Sutherland equation: $\tau_D=(R_C)^2/D$ <br>

*Note 3:* In simulation, there are effectively 2 time conversions: a conversion to real time, and an internal simulation timescale that can be calculated from the [length], [mass], and [energy] parameters; these value do NOT match because, due to numerical constraints, simulation uses an idealized relationship between momentum diffusivity (kinematic viscosity) and mass diffusivity. The resulting Schmidt number is on the order of 1, when in reality the Schmidt number is on the order of millions. <br>The take-away here is that in order to convert a constant that uses time into simulation units you should normalize it via energy scaling, and you should use the diffusion time to compare the timescale of colloidal simulations with experimental results.
