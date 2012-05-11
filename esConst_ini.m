nequil = 1000;
ntotal = 10000;

nangles = [2 3 4 5 6 7 8 9 10 15 20 30 50];
beta1s  = [1 30];
betaESs = [-10 -20 -30];
Vgss    = [0.6 1.5 2.5];

esConst_ini = cell(size(nangles,2), size(beta1s,2), size(betaESs,2), size(Vgss,2));

ic = 0;
for iangles = 1:size(nangles,2)
   for ibeta1 = 1:size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
            ic = ic + 1;
            disp(['starting run ',num2str(ic)]);
            C = TrajConfig;
            C.nangles        = nangles(iangles);
            C.Vgs            = Vgss(iVgs);
            C.beta1          = beta1s(ibeta1);
            C.betaES         = betaESs(ibetaES);
            energiesToSave   = 2;
            C.ESforces       = 0;
          
            C.tstep          = C.beta1;
            nsave1 = [1 1 1 0 0];
            t1 = TrajSegment(C, nsave1, energiesToSave, 0);
            t1.runTraj(ntotal);
            
            e1 = t1.ener(2,nequil+1:end);
            res.emean = mean(e1);
            figure(100);
            plot(e1);
            res.esd = 100*std(e1)/sqrt(ntotal-nequil);
            esConst_ini{iangles,ibeta1,ibetaES,iVgs} = res;
         end
      end
   end
end

%%

save('data_ES/esConst_ini.mat','esConst_ini');


