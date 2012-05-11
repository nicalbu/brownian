function res = temp(obj)
% returns temperature as a function of time

% k in kcal/mol/K =
%   1.3806504e-23 (J/mol K) * (1kJ/1000 J) * (1 cal/4.184 J) * 6.02e23 =
%
k = 0.00198649986;
nangles = obj.trajSegment(1).nangles;
nsteps = obj.nsteps('vels');
Inertia = obj.trajSegment(1).C.I;

% total kinetic energy versus time
res = zeros(1,nsteps);
vels = obj.vels;
for itime = 1:nsteps
   res(1,itime) = Inertia * ...
      (vels(:,itime)'*vels(:,itime))/nangles/k;
end
