# hippo_plos
README for Fortran Code

Overview

This Fortran code simulates the dynamics of cell growth, division, and interaction within a cellular environment. The code employs a computational model that integrates reaction kinetics, cellular mechanics, and stochastic processes to investigate the evolution of cellular systems over time. The model is parameterized and controlled by various constants and initial conditions defined in the parameter module.

Modules and Key Variables parameter Module The parameter module contains all the constants and parameters used throughout the simulation. Below are the key variables and their descriptions:

n_t, m_t, nlit, mlit: Define the total number of particles and timesteps. INTEGMX, INTEGST, INTEG: Parameters for random number generation. ldim, mdim: Dimensions of the cellular environment and mapping space. pi, cbrt2, dt: Mathematical constants and the time step used in the simulation. m0, k_fric, epsilon: Physical properties such as initial cell mass, frictional coefficient, and energy. al_LJ, lambda, gamma: Parameters for the Lennard-Jones potential and initial cell arrangement. Reaction constants: a1, a2, b1, a3 define the rates for the chemical reactions in the model. Xth, Xinit: Threshold concentration for cell cycle phases and initial concentrations. nrho, rhon, d_am: Parameters related to density points and color coding. Program HIPPO The main program HIPPO initializes and calls the subroutine sHIPPO, which handles the simulation.

Subroutine sHIPPO This is the core subroutine that:

Initializes cell properties (positions, velocities, sizes, etc.). Evolves the system over time using the Verlet integration method. Computes forces based on the Lennard-Jones potential. Handles cell growth, division, and reaction dynamics using the Runge-Kutta method. Outputs the simulation data at specified intervals. Subroutine sEVP Handles the Eigenvalue problem related to the reaction equations. It calculates the new concentration values based on the current state and updates them accordingly.

Utility Subroutines sAv: Multiplies a matrix by a vector. sinvA: Computes the inverse of a matrix. suniform_rn: Generates uniform random numbers using a linear congruential method. Compilation and Execution

Compilation: To compile the code, use a Fortran compiler (e.g., gfortran): sh> gfortran -o hippo_simulation your_code_file.f90
Execution: Run the compiled program: sh> ./hippo_simulation Input and Output Input The code does not require external input files. All initial conditions and parameters are hardcoded in the parameter module. Output The program generates several output files: surat_t231218_4.dat: Contains summary statistics like survival ratios over time. Xyzw_t231218_4.dat: Stores data related to cell positions, sizes, and other attributes. xyzXXXX.csv: Contains detailed data on the state of the system at specific time intervals, where XXXX is a timestep index. Customization To customize the simulation:
Modify the parameters in the parameter module to change initial conditions, time steps, or physical properties. Adjust the logic within sHIPPO for different growth models or interaction rules. License This code is intended for research and educational purposes. For other uses, please contact the author or institution.

Contact For any questions or further information, please reach out to the author at [umegaki@sigmath.es.osaka-u.ac.jp].
