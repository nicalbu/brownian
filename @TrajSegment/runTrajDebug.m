function res = runTrajDebug(obj,totalSteps)

for it = (obj.timeSteps+1):totalSteps
   obj.lastAngles = ones(obj.nangles,1) * it;
   obj.lastVels = ones(obj.nangles,1) * it;
   Egs = it;
   Eexc = ones(obj.nener-1,1) * it;
   obj.lastWf = ones(obj.nangles,obj.nwf) * it;
   cent = it;
   obj.store(it,obj.lastAngles,obj.lastVels,Egs,Eexc,obj.lastWf,cent);
end
obj.timeSteps = totalSteps;
