classdef TrajSegment < dynamicprops
   properties
      C                   % Configuration used for this segment
      % note that the reserved time in following arrays may be larger
      % than the actual time (actual time is timeSteps)
      angles              % (nangles, time)
      vels                % (nangles, time)
      nwf                 % number of wavefunctions to save
                          % nwf = -1 means save intensities in wf(:,1,time)
      wf                  % (nangles, nwf, time) excited state wfs
      cent                % (1, time)
      nener               % # ener to save; 1=GS, 2..  are excited states
      ener                % (nener,time)
      timeSteps=0         % number of time steps taken
      
      % For these, nsave = 10 means save 1,11,21,31 etc
      nsave               % (1,5) save data every nsave steps (0=nosave)
      % 1 = angles, 2 = vels, 3= ener, 4 = wf, 5 = cent
      % default is 1,0,0,0,0
      
      initialAngles       % (nangles,1) start trajectory with these angles
      %  Default is random values
      initialVels         % (nangles,1) start trajectory with these velocities
      %  Default is zeros
      initialAcc          % (nangles,1) used only for allen
      initialAcc2         % (nangles,1) used only for allen
      % These are saved to allow restart:
      lastAngles          % (nangles,1), angles
      lastVels             % (nangles,1), velocities (in radians/md units)
      lastAcc             % (nangles,1) to restart allen algorithm
      % to restart allen algoritm (acc from 2 steps ago)
      lastAcc2            % (nangles,1)
      lastWf              % (nangles,1), wavefunction
   end
   properties (Dependent, SetAccess = private)
      nangles;
   end
   methods (Static)
      function res = toIndex(type)
         type1 = validatestring(type,{'angles','vels','ener','wf','cent'});
         switch lower(type1)
            case 'angles'
               res = 1;
            case 'vels'
               res = 2;
            case 'ener'
               res = 3;
            case 'wf'
               res = 4;
            case 'cent'
               res = 5;
         end
      end
      function res = adjustTime(tInMD,units)
         timeUnits = validatestring(units,{'md','fs','ps'});
         switch lower(timeUnits)
            case ('fs')
               res = tInMD * 48.8;
            case ('ps')
               res = tInMD * 0.0488;
            case('md')
               res = tInMD;
         end
      end
   end
   methods (Static)
      [forcesGS, E] = forcesFromGS(V,angles,periodic)
      [forces,E,c,wf,flag] = forcesFromES(beta,angles,b,cen,width,cutoff)
   end
   methods (Access = private)
      res = runTrajDebug(obj,nsteps)
      res = runTrajAllen(obj,nsteps)
      res = runTrajLF(obj,nsteps)
      res = runTrajFricOnly(obj,nsteps)
   end
   methods
      function obj = TrajSegment(C,nsaveIn,nener,nwf)
         if (nargin < 1)
            C = TrajConfig;
         end
         if (nargin < 2)
            nsaveIn = [1 0 0 0 0]; % only save angles, and do so every step
         end
         if (nargin < 3)
            nener = 1;
         end
         if (nargin < 4)
            nwf = 1;
         end
         obj.C = C;
         obj.nsave = nsaveIn;
         obj.nener = nener;
         obj.nwf = nwf;
      end
      function res = get.nangles(obj)
         res = obj.C.nangles;
      end
      function setNsave(obj,value,type)
         if (obj.timeSteps > 0)
            error('cant change nsave once trajectory has started');
         end
         if (nargin < 3)
            type = 'angles';
         end
         index = TrajSegment.toIndex(type);
         obj.nsave(1,index) = value;
      end
      function res = getNsave(obj,type)
         if (nargin < 2)
            type = 'angles';
         end
         res = obj.nsave(1,TrajSegment.toIndex(type));
      end
      function res = time(obj,type,units)
         if (nargin < 2)
            type = 'angles';
         end
         if (nargin < 3)
            units = 'md';
         end
         index = TrajSegment.toIndex(type);
         nsave1 = obj.nsave(1,index);
         t1 = (0:nsave1:(obj.timeSteps-1)) * obj.C.tstep;
         res = TrajSegment.adjustTime(t1,units);
      end
      function res = totalTime(obj,units)
         if (nargin < 2)
            units = 'md';
         end
         t1 = (obj.timeSteps-1) * obj.C.tstep;
         res = TrajSegment.adjustTime(t1, units);
      end
      function res = data(obj,type)
         if (nargin < 2)
            type = 'angles';
         end
         index = TrajSegment.toIndex(type);
         size = obj.nsteps(type);
         switch (index)
            case 1
               res = obj.angles(:,1:size);
            case 2
               res = obj.vels(:,1:size);
            case 3
               res = obj.ener(:,1:size);
            case 4
               res = obj.wf(:,:,1:size);
            case 5
               res = obj.cent(:,1:size);
         end
      end
      function res = nsteps(obj,type)
         if (nargin < 2)
            type = 'angles';
         end
         index = TrajSegment.toIndex(type);
         res = floor((obj.timeSteps-1)/obj.nsave(1,index))+1;
         % not incredibly efficient, but definitely correct
         % nsave1 = obj.nsave(1,index);
         % res = size( (0:nsave1:(obj.timeSteps-1)), 2);
      end
      function reserve(obj,ntimeSteps)
         % ntimeSteps is total number of time steps (independent of nsave)
         % does nothing if requested space is < existing space
         if (obj.nsave(1,1) > 0) % angles
            needed = size( (0:obj.nsave(1,1):(ntimeSteps-1)), 2);
            if (needed > size(obj.angles,2))
               t1 = obj.angles;
               obj.angles = zeros(obj.nangles,needed);
               if (size(t1,2) > 0)
                  obj.angles(:,1:size(t1,2)) = t1;
               end
            end
         end
         if (obj.nsave(1,2) > 0) % vels
            needed = size( (0:obj.nsave(1,2):(ntimeSteps-1)), 2);
            if (needed > size(obj.vels,2))
               t1 = obj.vels;
               obj.vels = zeros(obj.nangles,needed);
               if (size(t1,2) > 0)
                  obj.vels(:,1:size(t1,2)) = t1;
               end
            end
         end
         if (obj.nsave(1,3) > 0) % ener
            needed = size( (0:obj.nsave(1,3):(ntimeSteps-1)), 2);
            if (needed > size(obj.ener,2))
               t1 = obj.ener;
               obj.ener = zeros(obj.nener,needed);
               if (size(t1,2) > 0)
                  obj.ener(:,1:size(t1,2)) = t1;
               end
            end
         end
         if (obj.nsave(1,4) > 0) % wf
            needed = size( (0:obj.nsave(1,4):(ntimeSteps-1)), 2);
            if (needed > size(obj.wf,3))
               t1 = obj.wf;
               obj.wf = zeros(obj.nangles,abs(obj.nwf),needed);
               if (size(t1,2) > 0)
                  obj.wf(:,:,1:size(t1,3)) = t1;
               end
            end
         end
         if (obj.nsave(1,5) > 0) % cent
            needed = size( (0:obj.nsave(1,5):(ntimeSteps-1)), 2);
            if (needed > size(obj.cent,2))
               t1 = obj.cent;
               obj.cent = zeros(2,needed);
               if (size(t1,2) > 0)
                  obj.cent(:,1:size(t1,2)) = t1;
               end
            end
         end
      end
      function extendFrom(obj,trajStart)
         obj.initialAngles = trajStart.lastAngles;
         obj.initialVels   = trajStart.lastVels;
         if (strcmp(obj.C.calcType, 'allen'))
            obj.initialAcc    = trajStart.lastAcc;
            obj.initialAcc2   = trajStart.lastAcc2;
         end
      end
      function res = runTraj(obj, totalSteps)
         % extend the trajectory to this # of total steps
         if (totalSteps > obj.timeSteps)
            % Reserve space (reserve does nothing if space already present)
            obj.reserve(totalSteps);
            % Initialize angles and velocites to random, or input values
            if (obj.timeSteps == 0)
               if (size(obj.initialAngles,1) == 0)
                  obj.lastAngles = (rand(obj.nangles,1)*pi/2);
               else
                  obj.lastAngles = obj.initialAngles;
               end
               if (size(obj.initialVels,1) == 0)
                  obj.lastVels = zeros(obj.nangles,1);
               else
                  obj.lastVels = obj.initialVels;
               end
               if (strcmp(obj.C.calcType,'allen'))
                  if (size(obj.initialAcc,1) == 0)
                     obj.lastAcc = zeros(obj.nangles,1);
                     obj.lastAcc2 = zeros(obj.nangles,1);
                  else
                     obj.lastAcc = obj.initialAcc;
                     obj.lastAcc2 = obj.initialAcc2;
                  end
               end
               % calculate energies and wf at this time step
               [~, Egs] = TrajSegment.forcesFromGS(obj.C.Vgs,obj.lastAngles);
               if (abs(obj.C.betaES) > 0)
                  [~,Eexc,c,wf1] = TrajSegment.forcesFromES(obj.C.betaES, ...
                     obj.lastAngles,[]);
                  if (obj.C.periodic)
                     nangles = obj.nangles;
                     sinAngles = sin( 2 * pi/nangles * (1:nangles) );
                     cosAngles = cos( 2 * pi/nangles * (1:nangles) );
                     sinAvg = sum(c.^2 .* sinAngles');
                     cosAvg = sum(c.^2 .* cosAngles');
                     cent1 = atan2(sinAvg,cosAvg);
                  else
                     cent1 = sum(c.^2 .* (1:obj.nangles)');
                  end
                  obj.lastWf  = c;
                  obj.store(1,obj.lastAngles, obj.lastVels, Egs,Eexc,wf1, cent1,0);
               else
                  wf1 = zeros(obj.nangles,obj.nangles);
                  obj.store(1,obj.lastAngles, obj.lastVels, Egs);
                  cent1 = 0;
               end
               obj.timeSteps = 1;
            end
            if (strcmp(obj.C.calcType,'debug'))
               obj.runTrajDebug(totalSteps);
            end
            if (strcmp(obj.C.calcType,'allen'))
               res = obj.runTrajAllen(totalSteps);
            end
         end
      end
   end % methods
end % class