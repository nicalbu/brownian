classdef Traj < dynamicprops
   properties
      trajSegment
   end
   methods
      function obj = Traj(trajSegment)
         if (nargin == 0)
            obj.trajSegment = TrajSegment.empty(0,0);
         else
            obj.trajSegment(1) = trajSegment;
         end
      end
      function res = nsegments(obj)
         % number of trajectory segments
         res = size(obj.trajSegment,2);
      end
      function addSegment(obj,trajSegment)
         obj.trajSegment(obj.nsegments + 1) = trajSegment;
      end
      function res = overlaps(obj,type)
         % res(i) = true if first data point of segment i
         %          overlaps the last data point of segment i-1
         if (nargin < 2)
            type = 'angles';
         end
         if (obj.nsegments > 0)
            res = zeros(obj.nsegments, 1);
            tstep = obj.trajSegment(1).C.tstep;
            for i=2:obj.nsegments
               if (obj.trajSegment(i-1).nsteps(type) > 0)
                  time1 = obj.trajSegment(i-1).time(type);
                  tlastdata = time1(1,size(time1,2));
                  tlast = obj.trajSegment(i-1).totalTime;
                  if (abs(tlast-tlastdata) < (1.0e-4*tstep) )
                     res(i) = true;
                  else
                     res(i) = false;
                  end
               end
            end
         else
            res = [];
         end
      end
      function res = nsteps(obj,type)
         if (nargin < 2)
            type = 'angles';
         end
         res = 0;
         olaps = obj.overlaps(type);
         for i=1:obj.nsegments
            if (olaps(i))
               res = res + obj.trajSegment(i).nsteps(type)-1;
            else
               res = res + obj.trajSegment(i).nsteps(type);
            end
         end
      end
      function res = time(obj,type,units)
         if (nargin < 2)
            type = 'angles';
         end
         if (nargin < 3)
            units = 'md';
         end
         res = zeros(1,obj.nsteps(type));
         olap = obj.overlaps(type);
         ic = 0;
         tlast = 0.0;
         for i = 1:obj.nsegments
            t1 = obj.trajSegment(i).time(type,units);
            if ( olap(i) )
               res(1,(ic+1):(ic+size(t1,2)-1)) = tlast ...
                  + t1(1,2:size(t1,2));
               ic = ic+size(t1,2)-1;
            else
               res(1,(ic+1):(ic+size(t1,2))) = tlast + t1;
               ic = ic + size(t1,2);
            end
            tlast = tlast + obj.trajSegment(i).totalTime(units);
         end
      end
      function res = data(obj,type)
         if (nargin < 2)
            type = 'angles';
         end
         if (obj.nsegments > 0)
            olap = obj.overlaps(type);
            t1 = obj.trajSegment(1).data(type);
            ndim = size(size(t1), 2);
            if (ndim == 2)
               res = zeros(size(t1,1),obj.nsteps(type));
            else
               res = zeros(size(t1,1),size(t1,2),obj.nsteps(type));
            end
            ic = 0;
            for i = 1:obj.nsegments
               t1 = obj.trajSegment(i).data(type);
               if ( olap(i) )
                  if (ndim == 2)
                     res(:,(ic+1):(ic+size(t1,2)-1)) = t1(:,2:size(t1,2));
                     ic = ic + size(t1,2)-1;
                  else
                     res(:,:,(ic+1):(ic+size(t1,3)-1)) = t1(:,:,2:size(t1,3));
                     ic = ic + size(t1,3)-1;
                  end
               else
                  if (ndim == 2)
                     res(:,(ic+1):(ic+size(t1,2))) = t1;
                     ic = ic + size(t1,2);
                  else
                     res(:,:,(ic+1):(ic+size(t1,3))) = t1;
                     ic = ic + size(t1,3);
                  end
               end
            end
         else
            res = [];
         end
      end
      function res = angles(obj)
         res = obj.data('angles');
      end
      function res = vels(obj)
         res = obj.data('vels');
      end
      function res = ener(obj)
         res = obj.data('ener');
      end
      function res = wf(obj)
         res = obj.data('wf');
      end
      function res = cent(obj)
         res = obj.data('cent');
      end
   end % methods
end