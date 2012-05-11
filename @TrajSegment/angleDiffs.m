function res = angleDiffs(obj)
% returns difference between angles at site i and i+1

nangles = obj.nangles;
nsteps = obj.nsteps('angles');
res = obj.angles(2:nangles,1:nsteps) - obj.angles(1:(nangles-1), 1:nsteps);
