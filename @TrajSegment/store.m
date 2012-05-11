function store(obj,it,angles,vels,Egs,Eexc,wf,cent,flag)

if ((obj.nsave(1,1) > 0) && (rem(it-1,obj.nsave(1,1)) == 0))
   ic = floor((it-1)/obj.nsave(1,1)) + 1;
   obj.angles(:,ic) = angles;
end
if ((obj.nsave(1,2) > 0) && (rem(it-1,obj.nsave(1,2)) == 0))
   ic = floor((it-1)/obj.nsave(1,2)) + 1;
   obj.vels(:,ic) = vels;
end
if ((obj.nsave(1,3) > 0) && (rem(it-1,obj.nsave(1,3)) == 0))
   ic = floor((it-1)/obj.nsave(1,3)) + 1;
   esave = zeros(obj.nener,1);
   esave(1,1) = Egs;
   if ((nargin > 4) && (obj.nener > 1))
      esave(2:obj.nener,1) = Eexc(1:(obj.nener-1),1);
   end
   obj.ener(:,ic) = esave;
end
if (nargin > 5)
   if ((obj.nsave(1,4) > 0) && (rem(it-1,obj.nsave(1,4)) == 0))
      ic = floor((it-1)/obj.nsave(1,4)) + 1;
      if (obj.nwf > 0)
         obj.wf(:,:,ic) = wf(:,1:obj.nwf);
      else
         % Assuming transition operator is equal on each unit cell, the
         % intensity is the sum of the wavefunction amplitudes. This is
         % done by summing over the first index, and saving the results in
         % the first column of obj.wf
         trans = sum(wf,1);
         obj.wf(:,1,ic) = trans(1,:);
      end
   end
end
if (nargin > 6)
   if ((obj.nsave(1,5) > 0) && (rem(it-1,obj.nsave(1,5)) == 0))
      ic = floor((it-1)/obj.nsave(1,5)) + 1;
      obj.cent(1,ic) = cent;
   end
end
if (nargin > 7)
   if ((obj.nsave(1,5) > 0) && (rem(it-1,obj.nsave(1,5)) == 0))
      ic = floor((it-1)/obj.nsave(1,5)) + 1;
      obj.cent(2,ic) = flag;
   end
end
