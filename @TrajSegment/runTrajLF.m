function res = runTrajLF(obj,totalSteps)
% Run trajectory 

% Simulation parameters
nangles    = obj.C.nangles;
tstep      = obj.C.tstep; 
gamma      = obj.C.beta1;
I          = obj.C.I;

% This is leap-frog integration from:
% W.F. Van Gunsteren and H.J.C. Berendsen, Molecular Simulation
% vol 1 pp 173-185.
% The process is spelled out on the bottom of page 180, and the steps 
% discussed below correspond to those listed there. 

%  kt -> kcal/mol
%  I  -> amu A^2
%  gamma -> md time units 

kT = 1.3806504*10^(-23)*temp*0.2390057*10^(-3)*6.02*10^23; % in kcal/mol
%    (J/K)(K) (kcal/J) * (mol)    = kT in kcal/mol

% XX at time t(n-1/2)
% to intialize, we sample from a gaussian distribution with zero mean and 
% width (3.18)
C323  = gamma * tstep - 3 + 4 * exp(-1.0*gamma*tstep/2.0) ...
   - exp(-1.0 * gamma *tstep);
eq318 = kT/(I*gamma^2) * C323;

assert(eq318 > 0.0, 'eq318 is negative');
width318 = sqrt(eq318);
XXplus = randn(nangles,1) * width318;

% Initialize variables used in the integration
D324  = 2 - exp(gamma * tstep/2.0) - exp(-1.0*gamma *tstep/2.0);
B325 = gamma * tstep * (exp(gamma*tstep)-1.0)...
   - 4.0 * (exp(gamma*tstep/2.0) -1.0)^2;
eq321 = (kT/I) * B325/C323;
assert(eq321 > 0.0, 'eq321 is negative');
width321 = sqrt(eq321);

eq329 = (kT/I)*(1.0 - exp(-gamma*tstep));
assert(eq329 > 0.0, 'eq329 is negative');
width329 = sqrt(eq329);

B325minus = -1.0*gamma * tstep * (exp(-1.0*gamma*tstep)-1.0)...
   - 4.0 * (exp(-1.0*gamma*tstep/2.0) -1.0)^2;
eq332 = -1.0 * (kT/I/gamma^2)*B325minus / ...
   (exp(-1.0*gamma*tstep)-1);
assert(eq332 > 0.0, 'eq332 is negative');
width332 = sqrt(eq332);

for istep = (obj.timeSteps+1):totalSteps
   % Step 2: Force evaluation
   if (C.type == 0)
      forces = 0.;
   else
      if (C.type == 2) && (istep > C.n1)
         V = C.V2;
      end
      [forcesGS,Egs]         = forcesFromGS(V,angles);
      forces                 = forcesGS;
      if (C.type == 3)
         [forcesES,Ees,newWF] = forcesFromCharge(beta,angles,ESwf);
%         figure(600);
%         [maxwf,imaxwf] = max(abs(newWF));
%         plot(newWF * sign(newWF(imaxwf)),'r-');
%         input "continue";
         xbar = sum( (newWF.^2).*(1:nangles)' ,1);
         x2mean = sum( (newWF.^2).*((1:nangles)').^2 ,1);
         wfsd = sqrt(x2mean - xbar^2);
         if (istep > C.n1)
            forces=forces+forcesES;
            ESwf = newWF;
         end
      end
   end
   % Step 3
   Yv =  randn(nangles,1) * width321;
   Vminus = XXplus * gamma * D324/C323 + Yv; % eq. 3.36
   Vplus  = randn(nangles,1) * width329;
   nextVel = vel * exp(-gamma * tstep) ...
      + (1/I)* forces*tstep * (1/gamma/tstep) * (1.0-exp(-gamma*tstep))...
      + Vplus - exp(-gamma*tstep) * Vminus;
   Yx = randn(nangles,1) * width332;
   % note that D324minus = D324
   XXminus = Vplus * (1/gamma) * D324/(exp(-gamma*tstep)-1) ...
      + Yx; % eq. 3.37
   XXplus = randn(nangles,1) * width318;
   % eq. 3.10
   nextAngles = angles + nextVel * (1.0/gamma) * ...
      (exp(gamma*tstep/2.0) - exp(-gamma*tstep/2.0)) ...
      + XXplus - XXminus;
   angles  = nextAngles;
   vel     = nextVel;
   
   if (rem(istep,nSave) == 0)
      isave               = isave + 1;
      anglesSave(:,isave) = angles(:);
      velSave(:,isave)    = vel(:);
      if (C.type == 3)
         EesSave(1,isave)    = Ees;
         wfsdSave(1,isave)   = wfsd;
      end
   end
end

res.angles = anglesSave;
res.vels   = velSave;
res.tstep  = tstep;
res.nSave  = nSave;
if (C.type == 3)
   res.EesSave = EesSave;
   res.wfsdSave = wfsdSave;
end
