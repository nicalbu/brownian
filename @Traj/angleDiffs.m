function res = angleDiffs(obj)
% returns difference between angles at site i and i+1

angles = obj.angles;
nangles = size(angles,1);
nsteps = size(angles,2);
res = angles(2:nangles,1:nsteps) - angles(1:(nangles-1), 1:nsteps);
