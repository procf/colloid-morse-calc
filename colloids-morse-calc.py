# colloids-morse-calc is a Python script for calculating the Morse Potential and related MD simulation parameters from experimental values

def morse_calc(c_d, M_d, c_Star, eta, rho_s, phi, d_c, rho_c, data_time, dt_Integration, period, sim_R_C, sim_eta0, sim_rho, sim_kT):
  # constants
  # ---------
  T = 298                                # room temperature [K]
  k_B = 1.381e-23                        # Boltzmann constant [J/K]
  kT = k_B*T                             # energy unit kT
  epsilon_0 = 8.854e-12                  # permittivity of free space [F/m = (C^2 s^2)/(kg m^3)]
  e_charge = 1.602e-19                   # electron charge [C]
  Na = 6.022e23                          # Avogadro's number [#/Mol]
  L_m3_conversion = 1000                 # 1000 L per m^3 (1L = 0.001m^3)
  g_real = 9.81                          # acceleration due to gravity (on Earth) [m/s^2]

  # calculate polymer radius of gyration from c_star
  # -----------------------------------------------
  r_g = (M_d/(Na*(4/3)*np.pi*c_star*1e3))**(1/3)           # polymer radius of gyration [m]
  print('DEPLETANT')
  print(' - depletant radius of gyration r_g =',round(r_g/1e-9,3),'nm')
  print(' - depletion attraction range 2*r_g =',round((2*r_g)/1e-9,3),'nm')

  print('\nCOLLOIDS')
  # calulate the colloid radius and deplant-colloid ratio
  # -----------------------------------------------------
  R_C = 0.5*d_c
  delta = r_g/R_C
  print(' - depletant:colloid ratio, delta =',round(delta,3))

  # calculate colloid mass (assuming perfect spherical particle)
  # ------------------------------------------------------------
  rho_c_m3 = rho_c * L_m3_conversion     # kg/L * 1000L/m^3 = [kg/m^3]
  V_C = (4/3)*np.pi*R_C**3               # volume of a perfectly spherical particle [m^3]
  mass_C = V_C * rho_c_m3                # estimated mass of one colloid [kg]
  print(' - approximate mass per colloid, mass_C =',round((mass_C*1000)/1e-9,4),'nanograms')

  # timescale
  # ---------
  # Use the Stokes–Einstein–Sutherland equation for diffusion of spherical particles 
  # through a liquid with low Reynolds number to approximate the colloid diffusion time
  # ---------
  D = kT / (6*math.pi*eta*R_C)           # Diffusion coefficient for colloids
  tau_D = round((R_C)**2 / D)            # bare particle diffusion time [s]

  print(' - approximate colloid diffusion time =',tau_D,'s')
  print(' - data collected at ~',round((data_time)/tau_D),'diffusion times')

  print('\nDEPLETION')
  # DEPLETION POTENTIAL
  # -------------------
  # Calcualte the Asakura-Oosawa-Vrij type Depletion Potential
  # for more information on this calculation, see:
  #   Lekkerkerker Ruinier Vis Colloids and the Depletion Interaction, 2nd Edition (https://doi.org/10.1007/978-3-031-52131-7) 
  #   Section 2.1 Depletion Interaction Due to Penetrable Hard Spheres, page 67
  # -------------------
  def depletion(r, R_C1, R_C2, r_g, c_d):
    
    n_d = c_d * L_m3_conversion * (1/M_d) * Na    # depletant bulk number density : [g/L] * [L/m3] * [Mol/g] * [#/Mol] = [#/m3]
    n_d_free = n_d/(1-phi)                        # depletant number density in available solvent volume
    Pi_dynamic = n_d_free #* kT                   # osmotic pressure [J/m^3] / kT = [kT/m^3]
    h = r - (R_C1 + R_C2)                         # surface-surface separation distance [m]
    
    overlap_distance = np.where(h < 0, 0, np.maximum(0, (2 * r_g) - h))  # zero overlap if h > 2r_g OR h < 0
    
    # standard form: R_i = R_j -> R_C ; (3 * (2*(R_C1*R_C2/(R_C1+R_C2))) + 2 * r_g + h / 2)
    # for polysdisperse systems: R_i != R_j -> R_C = <R_C> = 2*(R_C1*R_C2/(R_C1+R_C2)))
    overlap_volume = (np.pi / 6) * overlap_distance**2 * (3 * (2*(R_C1*R_C2/(R_C1+R_C2))) + 2 * r_g + h / 2)
    
    return -Pi_dynamic * overlap_volume                                   # [kT/m^3] * [m^3] = [kT]

  D0 = -depletion((2*R_C), R_C, R_C, r_g, c_d) #/ (k_B*T)
  print(' - attraction strength at particle contact (potential well minimum) : D0 =', round(D0,2),'kT')

  # gravity
  print('\nGRAVITY')
  # -------
  # Calculate a dimensionless gravitational Peclet number settling/resistance
  # -------
  rho_c_meters = rho_c * L_m3_conversion # kg/L -> kg/m^3
  rho_s_meters = rho_s * L_m3_conversion # kg/L -> kg/m^3
  delta_rho = rho_c_meters-rho_s_meters

  def GravPeclet(R_C, delta_rho, g_real, k_B, T):
    settling = (4/3)*np.pi*R_C**3*delta_rho*g_real
    resistance = (k_B*T)/R_C

    Pe = settling/resistance
    return(Pe)

  Pe_real = GravPeclet(R_C, delta_rho, g_real, k_B, T)
  F_G_real = Pe_real * (k_B*T/R_C)
  print(' - Gravitational Peclet number for colloids is',round(Pe_real,4))
  print(' - real F_G per colloid =',F_G_real,'N (kg m/s^2)')

  print('\n')
  # ----------------------
  # simulation conversions
  # ----------------------
  t1 = period * dt_Integration            # simulation time per frame [time]
  sim_V_C = (4/3)*np.pi*sim_R_C**3        # volume of a sim colloid 
  sim_mass_C = sim_V_C * sim_rho          # mass of a sim colloid [mass]

  # [length]
  sim_length = R_C                        # m/[length]

  # [energy]
  sim_energy = kT/sim_kT                  # J/[energy]

  # [mass]
  sim_mass = mass_C / sim_mass_C          # kg/[mass]

  # [time]
  sim_D = sim_kT/(6*math.pi*sim_eta0*sim_R_C) # diffusion coefficient = r^2 / tau
  sim_tau_D = (sim_R_C**2)/(sim_D)        # bare particle diffusion time [time]
  sim_time = tau_D / sim_tau_D            # s/[time]

  print('Conversions for [length], [energy], [mass], and [time] calculated\n')

  sim_time_tauD = t1/sim_tau_D                # diffusion time per frame
  s_per_frame = sim_time_tauD * tau_D         # s/[time] -- how many seconds per sim frame?
  #nframes = data_time / s_per_frame  # [time]   -- how many frames to match exp time?
  nframes = np.ceil(data_time/tau_D) / sim_time_tauD

  sim_R_C = round(R_C/sim_length)

  print('how much time per frame?           |',round(s_per_frame),'s =',round(sim_time_tauD,2),'diffusion times')
  print('how many frames to match exp time? |',int(np.ceil(nframes)),'frames =',round(np.ceil(nframes)*sim_time_tauD,2),'diffusion times')

  print('\n  - simulation colloid size, sim_R_C =',sim_R_C)

  # Define Morse Potential

  def morse(r, R_C1, R_C2, D_0, alpha):
    h_ij = r - (R_C1+R_C2)
    morse_potential = D_0 * (np.exp(-2 * alpha * h_ij) - 2 * np.exp(-alpha * h_ij))
    return morse_potential

  D0 = round(D0)
  print('  - D0 = U0/kT ~',D0)
  alpha = round(round((3/(r_g/sim_length))))
  print('  - alpha ~',alpha)


  # simulation gravity
  # ------------------
  # momentum works differently in simulation vs reality due to Schmidt number mis-match, so do NOT use the real time conversion
  # instead, use energy based scaling: J = kg * m^2 / s^2 --> m/s^2 = J/(kg*m) = [energy]/([mass]*[length])
  # NOTE: this is also how HOOMD-blue calculates "time": [time] = sqrt(([mass]*[length]**2)/[energy])
  g_sim = g_real / (sim_energy/(sim_mass*sim_length))
  print('  - g in sim units = ',round(g_sim,4))

  sim_delta_rho = delta_rho * (sim_length**3/sim_mass)

  def GravPeclet_sim(sim_R_C, sim_delta_rho, g_sim, sim_kT):
    settling = (4/3)*np.pi*sim_R_C**3*(sim_delta_rho)*g_sim
    resistance = sim_kT/sim_R_C

    Pe = settling/resistance
    return(Pe)

  Pe_sim = GravPeclet_sim(sim_R_C, sim_delta_rho, g_sim, sim_kT)
  print('  - Pe in sim units = ',round(Pe_sim,4),'(should be same as real Pe!)')
  F_G_sim = Pe_sim * sim_kT/1
  print('  - sim F_G per unit length =',round(F_G_sim,4))

  print('\nSCHMIDT NUMBER')
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

  print(' - real Schmidt number:',round(Sc_real_millions),'million')
  print(' - Schmidt number in simulation:',round(Sc_sim))
