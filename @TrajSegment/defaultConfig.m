function res = defaultConfig()
% Config sets all parameters for the simulation
% this routine returns a Config with default values

% variable that apply across calculations
res.nangles    = 200;
res.tstep      = 1.0;     % dynamics uses this time step (in md units)
res.nSave      = 1;       % save information every nSave time steps
res.temp       = 298;
res.beta1      = 1;
res.x          = 3.93586; %distance between thiophene unit, in angstroms
res.I          =91.1194;  %inertia per thiophene unit in amu*A^2

