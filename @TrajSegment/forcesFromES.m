function [forces,E,c,wf,flag] = forcesFromES(beta,angles,b,periodic,ESeps, ...
   cen,width,cutoff)
% INPUT: beta = coupling between adjacent rings when planar (should be <0)
%       angles = list of all angles of the rings
%       b = wavefunction of charge from previous time step
%           this is used to select the lowest energy wavefunction with
%           good overlap with the previous time step
%           if b==[] (i.e. size(b,1) == 0), then lowest state is taken 
%       peroidic = true for periodic boundary conditions (default: false)
%    The following two are optimization variables. We first diagonalize
%    the subblock around the cen position, and then see if charge is
%    localized here.
%         cen = expected center of charge
%         width = width to try first for finding charge
%         cutoff = if the norm of the wavefunction with 5 rings of the
%           boundary > cutoff, a full diagonalization is performed
% Output: F = force on each ring
%         E = excitation energies (lowest is the state used for forces)
%         c = excited state wavefunctions (ordered as in E)
%         flag = 1 if the optimization failed and 0 otherwise

if (nargin < 4)
   periodic = false;
end
if (nargin < 5)
    ESeps = 0.2;
end
% Using subblock optimization?
if (nargin < 6)
   optimization = false;
else
   optimization = true;
end

% from angles, we know the length of the chain (i.e # angles)
flag=0;
nangles= size(angles,1);
% Set up the Hamiltonian
% set up the entire hamiltonian
% Note that taking the absolute value of cos can be rationalized as
% choosing, as each ring is added, the phase of its orbital so that the
% overlap with the previous ring is positive. (May not work for periodic
% boundary conditions.)
hamil= zeros(nangles,nangles);
for i=1:nangles-1
   hamil(i,i+1) = beta*cos(angles(i+1)-angles(i));
   hamil(i+1,i) = hamil(i,i+1);
end
if (periodic)
   hamil(1,nangles) = beta*cos(angles(i+1)-angles(i));
   hamil(nangles,1) = hamil(1,nangles);
end

optimizationWorked = false;
if (optimization)
   % pull out the sub-block corresponding to range cen-width... cen+width
   %    but I need to make sure I don't go below 1 or above nangles
   cen = round(cen);
   lower_limit = max(1,cen-width);
   upper_limit = min(nangles,cen+width);
   % this is the range I want
   i1 = lower_limit:upper_limit;
   Hsub = hamil(i1,i1);
   [Vsub,D] = eig(Hsub);
   % I now have the wavefunctions in the range i1, stored in Vsub(i1,i1)
   % I'll copy this into the correct location of V(nangles, i1)
   % V will then hold the wavefunctions on the range 1..nangles
   % but I only have i1 states calculated so there are only i1 of these
   length_i1 = size(i1,2);
   V = zeros(nangles,length_i1);
   V(i1,:) = Vsub;
   
   if (size(b,1) > 0)
      % I want to make sure the charge hasn't jumped to some remote location
      % so I calculate the overlap of all of my wavefunctions with the wf
      % of the previous time step
      overlap = abs(V'*b);
      % If the lowest energy state has overlap > eps, I'll keep it
      % otherwise, I'll keep the lowest energy state with such overlap
      y = 1;
      if (max(overlap) < ESeps) % keep lowest if no state has big overlap
         y = 1;
      else
         while (overlap(y) < ESeps)
            y = y+1;
         end
      end
   else
      y = 1;
   end
   c=V(:,y);
   wf = V;
   E=diag(D);
   if (y>1)
      % need to swap 1 and y
      wf(:,1) = V(:,y);
      wf(:,y) = V(:,1);
      E(1) = D(y,y);
      E(y) = D(1,1);
   end
   
   % Test to see if the subblock approximation is ok
   % We'll look at the norm of the wavefunction within 5 of the boundaries
   % and compare to the cutoff parameter
   test = norm(c(lower_limit:lower_limit+5))+norm(c(upper_limit-5:upper_limit));
   if (test < cutoff)
      optimizationWorked = true;
   end
   
end % close of optimization

% if test fails or we aren't optimizing, do a full diagonalization
if (optimizationWorked == false)
   [V,D]=eig(hamil);
   % if b is zero, we will just take lowest energy state
   if (size(b,1) == 0)
      y = 1;
      flag = 1;
   else
      overlap = abs(V'*b);
%       if (overlap(1) > 0.5)
%          y=1;
%       else
%          [~,y]=max(overlap);
%       end
      y = 1;
      if (max(overlap) < ESeps) % keep lowest if no state has big overlap
         y = 1;
      else
         while (overlap(y) < ESeps)
            y = y+1;
         end
      end
      flag = y + overlap(y);
   end
   c=V(:,y);
   wf = V;
   E=diag(D);
   if (y>1)
      % need to swap 1 and y
      wf(:,1) = V(:,y);
      wf(:,y) = V(:,1);
      E(1) = D(y,y);
      E(y) = D(1,1);
   end

end

% Analytic derivatives of the charge energy
forces = zeros(nangles,1);
if (periodic)
   forces(1) = 2*beta*sin(angles(1)-angles(nangles)) *c(1)*c(nangles) ...
      -2*beta*sin(angles(2)-angles(1))*c(1)*c(2);
else
   forces(1) = -2*beta*sin(angles(2)-angles(1))*c(1)*c(2);
end
for i=2:nangles-1
   forces(i) = 2*beta*sin(angles(i)-angles(i-1))*c(i)*c(i-1) ...
      -2*beta*sin(angles(i+1)-angles(i))*c(i)*c(i+1);
end
if (periodic)
   forces(nangles) = 2*beta*sin(angles(nangles)-angles(nangles-1))...
      *c(nangles)*c(nangles-1) ...
      -2*beta*sin(angles(1)-angles(nangles))*c(1)*c(nangles);

else
   forces(nangles) = 2*beta*sin(angles(nangles)-angles(nangles-1))...
      *c(nangles)*c(nangles-1);
end