classdef TrajConfig
   % Configures the dynamics (potential energy, time step, brownian
   % parameters etc) 
   properties
      % Properties of the molecular system
      nangles    = 20      % # degrees of freedom in chain
      I          = 91.1194 % inertia per thiophene unit in amu*A^2
      periodic   = false   % periodic boundary conditions?
      % Potential energy configuration
      Vgs        = 1       % GS potential energy 
      betaES     = -20.0   % Excited state coupling between rings (kcal/mol)
      % Simulation parameters
      temp       = 298     % Temperature (in kelvin)    
      beta1      = 1       % Friction/random force parameter
      tstep      = 1.0     % dynamics uses this time step (in md units)
      calcType   = 'allen' % type of simulation (see set method)
      ESforces   = true    % include ES forces
      ESoverlap  = false   % if false, take lowest state
                           % if true, check overlap with previous time step
      ESeps      = 0.1     % threshold for overlap
      optWidth   = 0       % optimization parameter (no optimization if 0)
      optCutoff  = 0       % optimization parameter
   end
   methods (Static)
      function res = GSonly(nangles, Vgs)
         res = TrajConfig;
         res.Vgs = Vgs;
         res.nangles = nangles;
         res.betaES = 0.0;
      end
   end
   methods
      function obj = set.calcType(obj,type)
         obj.calcType = validatestring(type,{'allen','LF','fricOnly', ...
            'debug'});
      end
   end % methods
end % class