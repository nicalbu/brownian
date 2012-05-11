function [cfgs, ncfgs] = contents(obj)
% cfgs{:,1} = list of all unique configurations
% ncfgs{:,1} = array (n,1) which holds list of all file numbers
%              corresponding to this cfg

nfiles = obj.setNfile();

if (nfiles == 0)
   cfgs = [];
   ncfgs = [];
else
   % To speed things up, we'll read in all the Cfiles right at the start
   Cfiles = cell(nfiles,1);
   for ifile = 1:nfiles
      fileprefix = [obj.dataPath,'\',obj.filePrefix,'_',int2str(ifile)];
      load([fileprefix,'_cfg.mat'], 'Cfile');
      Cfiles{ifile} = Cfile;
   end
   files = 1:nfiles;
   icfg = 0;
   while (size(files,2) > 0)
      % take the first Cfile, and find all matches to this
      Ctarget = Cfiles{files(1)};
      res = [1];
      for i = 1:size(files,2)
         ifile = files(i);
         if (size(Library.comp_struct(Cfiles{ifile},Ctarget),1) == 0)
             res = [res; ifile];
         end
      end
      icfg = icfg+1;
      cfgs{icfg,1} = Ctarget;
      ncfgs{icfg,1} = res;
      files = setdiff(files,res'); % remove these files from the file list
      disp(['files left ', num2str(size(files,2))]);
   end
end