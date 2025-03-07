# colloids-morse-calc is a Python script for calculating the Morse Potential and related MD simulation parameters from experimental values

# - morse_calc
# - plot_morse_vs_depletion

import numpy as np
import matplotlib.pyplot as plt

# DEPLETION POTENTIAL
# -------------------
# Calcualte the Asakura-Oosawa-Vrij type Depletion Potential
# for more information on this calculation, see:
#   Lekkerkerker, Ruinier, Vis : Colloids and the Depletion Interaction, 2nd Edition (https://doi.org/10.1007/978-3-031-52131-7) 
#   Section 2.1 Depletion Interaction Due to Penetrable Hard Spheres, page 67
# -------------------
def depletion(r, R_C1, R_C2, r_g, c_d, M_d, phi):

  # constants
  L_m3_conversion = 1000   # 1000 L per m^3 (1L = 0.001m^3)
  Na = 6.022e23            # Avogadro's number [#/Mol]

  n_d = c_d * L_m3_conversion * (1/M_d) * Na    # depletant bulk number density : [g/L] * [L/m3] * [Mol/g] * [#/Mol] = [#/m3]
  n_d_free = n_d/(1-phi)                        # depletant number density in available solvent volume
  Pi_dynamic = n_d_free #* kT                   # osmotic pressure [J/m^3] / kT = [kT/m^3]
  h = r - (R_C1 + R_C2)                         # surface-surface separation distance [m]
    
  overlap_distance = np.where(h < 0, 0, np.maximum(0, (2 * r_g) - h))  # zero overlap if h > 2r_g OR h < 0
    
  # standard form: R_i = R_j -> R_C ; (3 * (2*(R_C1*R_C2/(R_C1+R_C2))) + 2 * r_g + h / 2)
  # for polysdisperse systems: R_i != R_j -> R_C = <R_C> = 2*(R_C1*R_C2/(R_C1+R_C2)))
  overlap_volume = (np.pi / 6) * overlap_distance**2 * (3 * (2*(R_C1*R_C2/(R_C1+R_C2))) + 2 * r_g + h / 2)
    
  return -Pi_dynamic * overlap_volume                                   # [kT/m^3] * [m^3] = [kT]


# MORSE POTENTIAL
# ---------------
# Define Morse Potential
# ---------------
def morse(r, R_C1, R_C2, D_0, alpha):
  h_ij = r - (R_C1+R_C2)
  morse_potential = D_0 * (np.exp(-2 * alpha * h_ij) - 2 * np.exp(-alpha * h_ij))
  return morse_potential


# GRAVITATIONAL PECLET NUMBER
# ---------------------------
# Dimensionless number for measuring the effect of gravity as compared to Brownian motion / gel attraction
# for more information on this calculation see:
#   Padmanabhan and Zia (2018) Soft Matter (https://doi.org/10.1039/C8SM00002F)
#   Russel, Saville, Schowalter : Colloidal Dispersions. Cambridge University Press (1989)
# ---------------------------
def GravPeclet(R_C, delta_rho, g_real, kT):
  settling = (4/3)*np.pi*R_C**3*delta_rho*g_real
  resistance = (kT)/R_C

  Pe = settling/resistance
  return(Pe)




###########
"""""""""""
MORSE_CALC
"""""""""""
###########

def morse_calc(c_d, M_d, c_star, eta, rho_s, phi, d_c, rho_c, data_time, dt_Integration, period, sim_eta0, sim_rho, sim_kT):

  # constants
  # ---------
  T = 298                  # room temperature [K]
  k_B = 1.381e-23          # Boltzmann constant [J/K]
  kT = k_B*T               # energy unit kT
  epsilon_0 = 8.854e-12    # permittivity of free space [F/m = (C^2 s^2)/(kg m^3)]
  e_charge = 1.602e-19     # electron charge [C]
  Na = 6.022e23            # Avogadro's number [#/Mol]
  L_m3_conversion = 1000   # 1000 L per m^3 (1L = 0.001m^3)
  g_real = 9.81            # acceleration due to gravity (on Earth) [m/s^2]


  #################
  ## REAL VALUES ##
  #################

  # calculate polymer radius of gyration from c_star
  # -----------------------------------------------
  r_g = (M_d/(Na*(4/3)*np.pi*c_star*1e3))**(1/3)  # [m]

  # calulate the colloid radius and deplant-colloid ratio
  # -----------------------------------------------------
  R_C = 0.5*d_c    # [m]
  delta = r_g/R_C  # dimensionless

  # calculate colloid mass (assuming perfect spherical particle)
  # ------------------------------------------------------------
  rho_c_m3 = rho_c * L_m3_conversion     # kg/L * 1000L/m^3 = [kg/m^3]
  V_C = (4/3)*np.pi*R_C**3               # volume of a perfectly spherical particle [m^3]
  mass_C = V_C * rho_c_m3                # estimated mass of one colloid [kg]

  # timescale
  # ---------
  # Use the Stokes–Einstein–Sutherland equation for diffusion of spherical particles 
  # through a liquid with low Reynolds number to approximate the colloid diffusion time
  # ---------
  D = kT / (6*np.pi*eta*R_C)           # Diffusion coefficient for colloids
  tau_D = round((R_C)**2 / D)            # bare particle diffusion time [s]

  # depletion attraction strength 
  # -----------------------------
  # the potential at colloid-colloid contact, AKA potential well minimum
  # -----------------------------
  D0 = -depletion((2*R_C), R_C, R_C, r_g, c_d, M_d, phi) #/ (k_B*T)

  # print real values
  print('DEPLETANT')
  print(' - depletant radius of gyration                | r_g =',round(r_g/1e-9,3),'nm')
  print(' - depletion attraction range                  | 2*r_g =',round((2*r_g)/1e-9,3),'nm')
  print(' - attraction strength (potential minimum)     | D0 =', round(D0,2),'kT')

  print('\nCOLLOIDS')
  print(' - depletant:colloid ratio                     | delta =',round(delta,3))
  print(' - approximate mass per colloid                | mass_C =',round((mass_C*1000)/1e-9,4),'nanograms')
  print(' - approximate colloid diffusion               | tau_D =',tau_D,'s')
  print(' - data collected at                           | t ~',round((data_time)/tau_D),'diffusion times')


  # -------
  # Calculate a dimensionless gravitational Peclet number settling/resistance
  # -------
  rho_c_meters = rho_c * L_m3_conversion # kg/L -> kg/m^3
  rho_s_meters = rho_s * L_m3_conversion # kg/L -> kg/m^3
  delta_rho = rho_c_meters-rho_s_meters
  Pe_real = GravPeclet(R_C, delta_rho, g_real, k_B*T)
  F_G_real = Pe_real * (k_B*T/R_C)

  print('\nGRAVITY')
  print(' - Gravitational Peclet number for colloids    | Pe_G = ',round(Pe_real,4))
  print(' - gravitational force on each colloid         | F_GC =',F_G_real,'N (kg m/s^2)')


  ################
  ## SIM VALUES ##
  ################

  # [length]
  sim_length = R_C                        # m/[length]
  sim_R_C = round(R_C/sim_length)         # colloid size in simulation units
  sim_V_C = (4/3)*np.pi*sim_R_C**3        # volume of a sim colloid 

  # [energy]
  sim_energy = kT/sim_kT                  # J/[energy]

  # [mass]
  sim_mass_C = sim_V_C * sim_rho          # mass of a sim colloid [mass]
  sim_mass = mass_C / sim_mass_C          # kg/[mass]

  # [time]
  t1 = period * dt_Integration            # simulation time per frame [time]
  sim_D = sim_kT/(6*np.pi*sim_eta0*sim_R_C)          # diffusion coefficient = r^2 / tau
  sim_tau_D = (sim_R_C**2)/(sim_D)                   # bare particle diffusion time [time]
  sim_time = tau_D / sim_tau_D                       # s/[time]
  sim_time_tauD = t1/sim_tau_D                       # diffusion time per frame
  s_per_frame = sim_time_tauD * tau_D                # s/[time] -- how many seconds per sim frame?
  nframes = np.ceil(data_time/tau_D) / sim_time_tauD # [time]   -- how many frames to match exp time?

  # calculate Morse Potential
  # -------------------------
  D0_M = round(D0)
  alpha = round(round((3/(r_g/sim_length))))

  print('\nMORSE POTENTIAL')
  print(' attraction strength                           | D0 = U0/kT ~',D0_M)
  print(' attraction range parameter                    | alpha ~',alpha)

  # simulation gravity
  # ------------------
  # momentum works differently in simulation vs reality due to Schmidt number mis-match, so do NOT use the real time conversion
  # instead, use energy based scaling: J = kg * m^2 / s^2 --> m/s^2 = J/(kg*m) = [energy]/([mass]*[length])
  # NOTE: this is also how HOOMD-blue calculates "time": [time] = sqrt(([mass]*[length]**2)/[energy])
  g_sim = g_real / (sim_energy/(sim_mass*sim_length))
  sim_delta_rho = delta_rho * (sim_length**3/sim_mass)
  Pe_sim = GravPeclet(sim_R_C, sim_delta_rho, g_sim, sim_kT)
  F_G_sim = Pe_sim * sim_kT/1

  # print additional simulation parameters
  print('\nSIMULATION PARAMETERS')
  print('  how much time per frame?                     | 1 frame =',round(s_per_frame),'s =',round(sim_time_tauD,2),'diffusion times')
  print('  how many frames to match exp time?           | exp recorded at',int(np.ceil(nframes)),'frames =',round(np.ceil(nframes)*sim_time_tauD,2),'diffusion times')
  print('  - simulation colloid size                    | R_C =',sim_R_C)
  print('  - timestep                                   | dt_Integration =',dt_Integration)
  print('  - sim recording interval                     | period =',period)
  print('  - viscosity parameter                        | eta0 =',sim_eta0)
  print('  - system number density                      | rho =',sim_rho)
  print('  - system temperature                         | kT=',sim_kT)
  print('  - gravitational acceleration                 | g_sim =',round(g_sim,4))
  print('  - gravitational Peclet number                | Pe_G =',round(Pe_sim,4),'(should be same as real Pe_G!)')
  print('  - F_G per unit length                        | F_G =',round(F_G_sim,4))

  print('\nSIMULATION UNITS:')
  print('  - simulation unit length                     |',sim_length,'[m]')
  print('  - simulation unit mass                       |',sim_mass,'[kg]')
  print('  - simulation unit energy                     |',sim_energy,'[J]')
  print('  - simulation unit time (from diffusion time) | ',sim_time,'[s]')


  print('\nSCHMIDT NUMBER COMPARISON')
  # Schmidt number
  # --------------
  # for numerical efficienty, momentum is not the same in simulation as it is in experiment;
  # therefore, the Schmidt number (the ratio of momentum diffusivity over mass diffusivity) is always
  # significantly different than in experiment
  # You can check this for yourself:
  # Schmidt Number = (kinematic viscosity $\nu$ / Diffusion coefficient D) = dynamic visocisty $\mu$ / (fluid density $\rho_s$ * D)
  Sc_real = eta / (rho_s_meters * D)
  Sc_real_millions = Sc_real / 1e6

  Sc_sim = sim_eta0 / ((rho_s_meters*(sim_length**3/sim_mass))*sim_D)

  print(' - real Schmidt number                         |',round(Sc_real_millions),'million')
  print(' - Schmidt number in simulation                |',round(Sc_sim))

  print('\n')

  return(R_C, r_g, D0_M, alpha)



#######################
"""""""""""""""""""""""
PLOT_MORSE_V_DEPLETION
"""""""""""""""""""""""
#######################


def plot_morse_v_depletion(R_C, r_g, c_d, M_d, phi, D_0, alpha):
  # PLOT DEPLETION AND MORSE POTENTIAL
  r_min = 1e-11
  r_max = 1.2e-6 # plot range (must extend beyond R_C_bare+R_C_bare+2*r_g)
  r = np.linspace(r_min, r_max, 100_000)
  bare_contact = 2*R_C
  #colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

  # Potential calculations
  D_0 = -depletion(bare_contact, R_C, R_C, r_g, c_d, M_d, phi)
  U_morse = morse(r, R_C, R_C, D_0, alpha/R_C)
  U_depl = depletion(r, R_C, R_C, r_g, c_d, M_d, phi)

  ### PLOT MORSE vs. DEPLETION POTENTIAL
  plt.figure(figsize=(8, 6))

  plt.title(r"Depletion vs. Morse", fontsize=18)
  plt.xlabel("$h_{ij}/R_C$", fontsize=18)
  plt.ylabel(r"$U / k_B T$", fontsize=18)

  # add reference lines
  plt.axhline(0, color='black', linewidth=0.8, linestyle='--') # zero interaction
  plt.axvline(0, color='black', linewidth=0.8, linestyle='--') # colloid-colloid contact

  # label colloid zone
  plt.axvspan(-r_max / 1e-9, 0, alpha=0.05, color='grey')
  plt.text(0.08, 0.05, '$colloid$', fontsize=12,
              horizontalalignment='left', verticalalignment='bottom', 
              transform=plt.gca().transAxes)

  # include parameter label
  experimental_param = f"{c_d}mg/mL, $r_g=${round(r_g/1e-9,2)}nm"
  morse_param = f"$U_0/kT={round(D_0)}, \\alpha={round(alpha)}$"
  #props = dict(boxstyle='square', facecolor='grey', alpha=0.05)
  #plt.text(0.95, 0.2, f"{experimental_param}\n{morse_param}", transform=plt.gca().transAxes, 
  #       fontsize=14, verticalalignment='top', horizontalalignment='right', multialignment='center', bbox=props)
  plt.text(0.95, 0.15, experimental_param, transform=plt.gca().transAxes, 
         fontsize=14, verticalalignment='top', horizontalalignment='right', 
         color='black')
  plt.text(0.95, 0.09, morse_param, transform=plt.gca().transAxes, 
         fontsize=14, verticalalignment='top', horizontalalignment='right',
         color='#e66c00')

  plt.plot((r-bare_contact) / R_C, U_depl, ls='-', lw='1.5', color='black', label='Depletion')
  plt.plot((r-bare_contact) / R_C, U_morse, ls='-', lw='1.5', color='#ff7800', label='Morse') 

  #plt.legend(fontsize=14, loc='center left', bbox_to_anchor=(1,0.5))
  plt.legend(fontsize=14, loc='upper right')

  plt.xlim(-0.05, 0.75*(r_max-bare_contact)/R_C)
  plt.ylim((-1.25*D_0),(0.5*D_0))

  #plt.grid(True)

  plt.show()

