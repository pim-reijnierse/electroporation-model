%% MODEL PARAMETERS

% Time stepping
par.tlist1a = '10^{range(-8,(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))} tpulse+10^{range(-8,(log10(1e-4)+8)/50,log10(1e-4))} tpulse+10^{range(log10(1e-4),(log10(2e2)-log10(1e-4))/350, log10(2e2))}';
par.tlist1b = '10^{range(-8,(log10(tpulse[1/s])+8)/20,log10(tpulse[1/s]))}  tpulse+10^{range(-8,(log10(1e-4)+8)/5,log10(1e-4))}  tpulse+10^{range(log10(1e-4),(log10(2e2)-log10(1e-4))/50, log10(2e2))}';

% Pulse parameters and temperature
par.Eapp = 1e5; % Applied electric field [V/m]
par.tpulse = 1e-3; % Pulse duration [s]
par.trise = 1e-6; % Pulse rise/fall time [s]
par.npulse = 1; % Number of pulses
par.fp = 1e-10; % Repetition frequency [Hz]
par.T = 293; % Temperature [K]

% Geometrical parameters
par.rcell = 8.55e-6; % Cell radius [m]
par.dcell = 2*par.rcell; % Cell height [m] (for cylindrical cells)
par.del = 20*par.rcell; % Radius of the extracellular domain [m] 

% Cell electrical parameters
par.sigma_e = 1.58; % External conductivity [S/m]
par.sigma_i = 0.3; % Internal conductivity [S/m]
par.epsilon_e = 72; % External permittivity []
par.epsilon_i = par.epsilon_e; % Internal permittivity []
par.eps0 = 8.8542e-12; % Vacuum permittivity [F/m]
par.dm = 5e-9; % Membrane thickness [m]
par.Cm = par.eps0*5 / par.dm; % Membrane capacitance [F/m^2]
par.Gm = 9.5e-9 / par.dm; % Membrane conductance [S/m^2]
par.Urest = -50e-3; % Resting voltage (Um = Vi - Ve)

% Parameters for molecular transport
% Idx 1 = Lucifer Yellow
par.Nsolutes = 1; % Number of solute molecules simulated; Choose between 1 or 3
par.r0 = 0.19e-9; % Radius of charge carrier
par.l0 = 2*par.r0; % Length of charge carrier
par.r1 = 0.61e-9; % Radius of solute 1
par.l1 = 1.46e-9; % Length of solute 1
par.z1 = -2; % Valence of solute 1
par.D1_e = 4.77e-10; % Extracellular diff. coeficient of solute 1 [m^2/s]
par.D1_i = 0.25*4.77e-10; % Intracellular diff. coeficient of solute 1 [m^2/s]
par.c1_i0 = 0; % Initial intracellular concnetration of solute 1 [mol/m^3]
par.c1_e0 = 1; % Initial extracellular concnetration of solute 1 [mol/m^3]

par.r2 = 0.61e-9; % Radius of solute 2
par.l2 = 1.46e-9; % Length of solute 2
par.z2 = 1; % Valence of solute 2
par.D2_i = 0; % Intracellular diff. coeficient of solute 2 [m^2/s]
par.D2_e = 0; % Extracellular diff. coeficient of solute 2 [m^2/s]
par.c2_i0 = 0; % Initial intracellular concnetration of solute 2 [mol/m^3]
par.c2_e0 = 0; % Initial extracellular concnetration of solute 2 [mol/m^3]
par.r12 = 0.61e-9; % Radius of bound solutes 1&2
par.l12 = 1.46e-9; % Length of bound solutes 1&2
par.z12 = 1; % Valence of bound solutes 1&2
par.D12_i = 0; % Intracellular diff. coeficient of bound solutes 1&2 [m^2/s]
par.D12_e = 0; % Extracellular diff. coeficient of bound solutes 1&2 [m^2/s]
par.c12_i0 = 0; % Initial intracellular concnetration of bound solutes 1&2 [mol/m^3]
par.c12_e0 = 0; % Initial extracellular concnetration of bound solutes 1&2 [mol/m^3]
par.kass = 0; % Association constant [1/(s*mol/m^3)]
par.kdis = 0; % Dissociation constant [1/s]
% 1 M = 1 mol/l = 1e3 mol/m^3;

% Depth for convolved image
par.y_opt = 18e-6;