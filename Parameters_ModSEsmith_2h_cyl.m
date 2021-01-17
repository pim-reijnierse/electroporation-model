%% ELECTROPORATION PARAMETERS

% Electroporation parameters
par.rstar = 0.65e-9; % Pore radius at hydrophobic and hydrophilic transition [m]
par.rpmin = 1.0e-9; % Local minimum in the pore radius space at Um = 0 [m]
par.rpmax = 12e-9; % Maximum pore radius [m]
par.alpha = 1e9; % Pore creation rate density [1/(m^2*s)]
par.Uep1 = 1/11; % Characteristic voltage of electroporation (asymmetric) [V]
par.Uep2 = sqrt(1/18); % Characteristic voltage of electroporation (symmetric) [V]
par.Dp = 2e-13; % Diffusion coefficient in the pore radius space [m^2/s]
par.gamma = 2e-11; % Edge tension [N]
par.Gamma0 = 1e-5; % Surface tension of nonelectroporated membrane [N/m]
par.Gamma1 = 2e-2; % Hydrocarbon-water interfacial tension [N/m]
par.Fmax = 6.9e-10; % Parameter in the electrical force tending to expand the pore [N/V^2]
par.rh = 0.95e-9; % Parameter in the electrical force tending to expand the pore [m]
par.rt = 0.23e-9; % Parameter in the electrical force tending to expand the pore [m]
par.B = 1.6301e-19; % Parameter in Wster pore energy [J]
par.b = 3.5341; % Parameter in Wster pore energy []
par.nu = 0.25; % Relative length of the pore entrance
par.fprot = 0.5; % Areal fraction of proteins
par.taup = 4; % Pore resealing time constant [s]

% Options
par.expr_molflux = 'Smith'; % Choose between: 'Li', 'Smith'
par.expr_Wsurf = 'Smith';   % Choose between: 'Li', 'Smith', 'ConstantTension'
par.expr_Gp = 'Smith';      % Choose between: 'Li', 'Smith', 'Toroidal'

% Discretization
par.axi2D = false;
par.Ntheta = 60; % Number of discretization angles
par.dcomp2 = 5e-9; % Separation between lines in Component 2
