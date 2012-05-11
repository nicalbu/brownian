classdef Library < dynamicprops
   properties (SetAccess = private)
      dataPath     % directory for the data (without a trailing \)
      filePrefix   % prefix used fo cfg and data files
      nfile;       % current number of files (-1 if not net set)
   end
   methods (Static)
      [fs1, fs2, er] = comp_struct(s1,s2,n1,n2,p,tol);
   end
   methods
      function obj = Library(dataPathIn,filePrefixIn)
         if (nargin < 1)
            dataPathIn = 'data';
         end
         if (nargin < 2)
            filePrefixIn = 'lib';
         end
         obj.dataPath = dataPathIn;
         obj.filePrefix = filePrefixIn;
         obj.nfile = -1;
         % create the directory if it doesn't exist
         if (exist(obj.dataPath,'dir') ~= 7)
            status = mkdir(obj.dataPath);
            if (status == 0)
               error(['unable to create directory',obj.dataPath]);
            end
         end
      end
      function res = setNfile(obj)
         ifile = 0;
         failed = false;
         while (~failed)
            ifile = ifile + 1;
            fileprefix = [obj.dataPath,'\',obj.filePrefix,'_',int2str(ifile)];
            cfgfilename = [fileprefix,'_cfg.mat'];
            try % try to open the file
               load(cfgfilename,'Cfile');
            catch
               failed = true;
            end
         end
         obj.nfile = ifile -1;
         res = obj.nfile;
      end
      function res = store(obj,Cin,resFile)
         if (obj.nfile == -1)
            % figure out current number of files
            obj.setNfile();
         end
         obj.nfile = obj.nfile + 1;
         fileprefix = [obj.dataPath,'\',obj.filePrefix,'_',int2str(obj.nfile)];
         Cfile = Cin;
         save([fileprefix,'_cfg.mat'],  'Cfile' );
         save([fileprefix,'_calc.mat'], 'resFile' );
         res = obj.nfile;
      end
      function resFile = retrieve(obj, ifile)
         fileprefix = [obj.dataPath,'\',obj.filePrefix,'_',int2str(ifile)];
         load([fileprefix,'_calc.mat'], 'resFile' );
      end
      function res = find(obj, Ctarget, findAll)
         % returns file numbers that match Ctarget
         res = [];
         if (nargin < 3)
            findAll = true;
         end
         ifile = 0;
         nfound = 0;
         finished = false;
         failed = false;
         while (~finished)
            ifile = ifile + 1;
            fileprefix = [obj.dataPath,'\',obj.filePrefix,'_',int2str(ifile)];
            %disp(['looking for ',fileprefix]);
            cfgfilename = [fileprefix,'_cfg.mat'];
            try % try to open the file
               load(cfgfilename,'Cfile');
            catch
               failed = true;
               finished = true;
               %disp('file not found');
            end
            if (~failed)
               if (size(Library.comp_struct(Cfile,Ctarget),1) == 0)
                  %disp('file matches');
                  nfound = nfound+1;
                  if (findAll)
                     res(nfound,1) = ifile;
                  else
                     res = ifile;
                     finished = true;
                  end
               %else
               %   disp('file does not match');
               end
            end
         end
      end
   end
end
