nequil = 1000;
ntotal = 10000;

nangles = [2 3 4 5 6 7 8 9 10 15 20 30 50];
beta1s  = [1 30];
betaESs = [-10 -20 -30];
Vgss    = [0.6 1.5 2.5];

esConst_inf = cell(size(nangles,2), size(beta1s,2), size(betaESs,2), size(Vgss,2));

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
          
            C.tstep          = C.beta1;
            nsave1 = [1 1 1 0 0];
            t1 = TrajSegment(C, nsave1, energiesToSave, 0);
            t1.runTraj(ntotal);
            
            e1 = t1.ener(2,nequil+1:end);
            res.emean = mean(e1);
            figure(100);
            plot(e1);
            res.esd = 100*std(e1)/sqrt(ntotal-nequil);
            esConst_inf{iangles,ibeta1,ibetaES,iVgs} = res;
         end
      end
   end
end

%%

save('data_ES/esConst_inf.mat','esConst_inf');


%%

load('data/esConstant.mat');

%%
ic = 0;
for iangles = 1:size(nangles,2)
   for ibeta1 = 1:size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
             ic = ic + 1;
             a(ic) = esConstant{iangles,ibeta1,ibetaES,iVgs}.esd;
         end
      end
   end
end
             
             
%%
load('data/rawES.mat');
load('data/esConstant.mat');

%%

i1 = 13;
i2 = 2;
i3 = 4;
i4 = 3;

time = rawES{i1,i2,i3,i4}.time;
ener = rawES{i1,i2,i3,i4}.ener;

time0 = time(201:end) - time(201);
ener0 = ener(2,201:end);

y_h = esConstant{i1,i2,i3,i4}.emean;
r0 = time0(299);
y_test = zeros(1,size(time0,2));
y1     = zeros(1,size(time0,2));
y2     = zeros(1,size(time0,2));


for i = 1:size(time0,2)
    y1(i) = (r0/time(200+i))^12; %exp(-time(200+i)/r0);
    y2(i) = - 2*(r0/time(200+i))^6;
    y_test(i) =  -20.0*(y1(i) + y2(i));
end

%%

for iangles = 1:size(nangles,2)
   for ibeta1 = 1:size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
             time = rawES{iangles,ibeta1,ibetaES,iVgs}.time;
             ener = rawES{iangles,ibeta1,ibetaES,iVgs}.ener;
             
             time0 = time(201:end) - time(201);
             ener0 = ener(2,201:end);
             
             y_h = esConstant{iangles,ibeta1,ibetaES,iVgs}.emean;
             
             disp(['plotting nangles ',num2str(nangles(iangles)), ...
                    ' beta1 ', num2str(beta1s(ibeta1)), ...
                    ' betaES ', num2str(betaESs(ibetaES)), ...
                    ' Vgs ', num2str(Vgss(iVgs)), ...
                    ' y_h ', num2str(y_h)]);
             
             figure(100)
             %plot(time0,y_test,'b');
             hold on
             plot([time0(1);time0(end)],[y_h,y_h],'r');
             hold on
             plot(time0,ener0,'b');
             
             pause
         end
      end
   end
end


%%
figure(200)
plot(y1,'ro')
hold on
plot(y2,'bo')
hold on
plot(y_test,'go')
             