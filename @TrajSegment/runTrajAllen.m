function res = runTrajAllen(obj,totalSteps)
% Run trajectory (res is number of times optimization failed) 
res = 0;
% Simulation parameters
nangles    = obj.C.nangles;
tstep      = obj.C.tstep; 
beta1      = obj.C.beta1;
I          = obj.C.I;
periodic   = obj.C.periodic;
gamma = beta1 * tstep;
kT = 1.3806504*10^(-23) * obj.C.temp;

% Constants used in the simulation
c0    = exp(-gamma);
c1    = (1-c0)/gamma;
c2    = (1-c1)/gamma;
c3    = (1/2-c2)/gamma;
ca    = c1;
cb    = c2 + c3;
cc    = -c3;
cd    = c0;
ce    = c2 - c0 * c3/c1;
cf    = c1 - c2 + 2 * c0 * c3/c1;
cg    = -c0 * c3/c1;

muBr = 0;
muBv = 0;

% periodic  1..N  goes from 2Pi  angle = 2*Pi/N * (1:N)
if (periodic)
   sinAngles = sin( 2 * pi/nangles * (1:nangles) );
   cosAngles = cos( 2 * pi/nangles * (1:nangles) );
end

angles = obj.lastAngles;
vels = obj.lastVels;
acc = obj.lastAcc;
lastAcc  = obj.lastAcc2;
if (obj.C.betaES ~= 0)
   wf = obj.lastWf;
   if (periodic)
      sinAvg = sum(wf.^2 .* sinAngles');
      cosAvg = sum(wf.^2 .* cosAngles');
      cent = atan2(sinAvg,cosAvg);
   else
      cent = sum(wf.^2 .* (1:obj.nangles)');
   end
else
   wf = [];
   cent = 0.0;
end
for istep = (obj.timeSteps+1):totalSteps
   % Eq. A5 from Allen, Molecular Physics vol. 40, 1073-1087 (1980).
   
   if (kT ~= 0)
      sigmasquareBr = 2*beta1*kT*0.2390057*10^(-3)* ...
         6.02*10^23/I*(tstep/beta1^2-2/beta1^3*(1-exp(-beta1*tstep)) ...
         +1/2/beta1^3*(1-exp(-2*beta1*tstep)));
      sigmasquareBv = 2*beta1*kT*0.2390057*10^(-3)* ...
         6.02*10^23/I*(1/2/beta1*(1-exp(-2*beta1*tstep)));
      rho           = 2*beta1*kT*0.2390057*10^(-3)*...
         6.02*10^23/I*(1/beta1^2*(1-exp(-beta1*tstep)) ...
         -1/2/beta1^2*(1-exp(-2*beta1*tstep)));
      mu            = [muBr,muBv];
      SIGMA         = [sigmasquareBr,rho;rho,sigmasquareBv];
      B             = mvnrnd(mu,SIGMA,nangles);
      Br            = B(:,1);
      Bv            = B(:,2);
   else
      Br            = zeros(nangles,1);
      Bv            = zeros(nangles,1);
   end
      
   % Eq. A3 from Allen, Molecular Physics vol. 40, 1073-1087 (1980).
   nextAngles = angles + ca * tstep * vels + cb * tstep^2 * acc + ...
      cc * tstep^2 * lastAcc + Br;
   
   % force evaluation
   [forces,Egs] = TrajSegment.forcesFromGS(obj.C.Vgs, nextAngles, ...
      periodic);
   
   if (obj.C.betaES ~= 0) % need to do ES calc
      if (obj.C.ESoverlap == false) % set wf to [] to avoid overlap calc
         wf = [];
      end
      if (obj.C.optWidth > 0)
         [forcesES,Eexc,wf,ESwf,flag] = ...
            TrajSegment.forcesFromES(obj.C.betaES, nextAngles, wf, ...
            periodic, obj.C.ESeps, cent,obj.C.optWidth,obj.C.optCutoff);
         if (flag)
            res = res + 1;
         end
      else
         [forcesES,Eexc,wf,ESwf,flag] = ...
            TrajSegment.forcesFromES(obj.C.betaES, nextAngles, wf,...
            periodic, obj.C.ESeps);
        if (flag)
            res = res + 1;
        end
      end
      if (periodic)
         sinAvg = sum(wf.^2 .* sinAngles');
         cosAvg = sum(wf.^2 .* cosAngles');
         cent = atan2(sinAvg,cosAvg);
      else
         cent = sum(wf.^2 .* (1:obj.nangles)');
      end
      if (obj.C.ESforces)
         forces = forces + forcesES;
      end
   end
   
   nextAcc = forces/I;
   
   % Eq. A4 from Allen, Molecular Physics vol. 40, 1073-1087 (1980).
   nextVel = cd * vels + ce * tstep * nextAcc + ... 
      cf * tstep * acc + cg * tstep * lastAcc + Bv;
   
   angles  = nextAngles;
   vels     = nextVel;
   lastAcc = acc;
   acc     = nextAcc;
   
   if (abs(obj.C.betaES) > 0.0)
      obj.store(istep,angles,vels,Egs,Eexc,ESwf,cent,flag);
   else
      obj.store(istep,angles,vels,Egs);
   end
   
end

obj.lastAngles = angles;
obj.lastVels = vels;
obj.lastAcc = acc;
obj.lastAcc2 = lastAcc;
obj.lastWf = wf;
obj.timeSteps = totalSteps;

