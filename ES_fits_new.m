% fitting raw ES data with new method
% in which c is constrained for t --> Inf,
% and initial amplitudes are also constrained.

clear classes;

nangles = [2 3 4 5 6 7 8 9 10 15 20 30 50];
beta1s  = [1 30];
betaESs = [-10 -20 -30];
Vgss    = [0.6 1.5 2.5];

rawES = cell(size(nangles,2),size(beta1s,2),size(betaESs,2),size(Vgss,2),5);
ESres = cell(size(nangles,2),size(beta1s,2),size(betaESs,2),size(Vgss,2),2); 

Clib.type = 'ediffuse';
Clib.nangles = 2;
Clib.Vgs = 0.6;
Clib.betaES = -10;
Clib.beta1 = 1;
%Nsave: 1 = angles, 2 = vels, 3= ener, 4 = wf, 5 = cent
Clib.time1 = 1; % time step will be beta1 * Clib.time1
Clib.nsave1  = [1 1 1 0 0];
Clib.nsteps1 = 200;

Clib.time2 = 0.1;
Clib.nsave2  = [1 1 1 0 0];
Clib.nsteps2 = 1501;

%%

load('data_ES/rawES_final.mat');
load('data_ES/esConst_inf.mat');
load('data_ES/esConst_ini.mat');

%% New Fit for ES

for iangles = 1:size(nangles,2)
   for ibeta1 = 1:size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
             Clib.nangles = nangles(iangles);
             Clib.beta1   = beta1s(ibeta1);
             Clib.betaES  = betaESs(ibetaES);
             Clib.Vgs     = Vgss(iVgs);
             
             disp(['loading nangles ',num2str(Clib.nangles), ...
                 ' beta1 ', num2str(Clib.beta1), ...
                 ' betaES ', num2str(Clib.betaES), ...
                 ' Vgs ', num2str(Clib.Vgs)] );
             cESinf = esConst_inf{iangles,ibeta1,ibetaES,iVgs}.emean;
             cESini = esConst_ini{iangles,ibeta1,ibetaES,iVgs}.emean;
             araw   = rawES{iangles,ibeta1,ibetaES,iVgs};
             time       = araw.time;
             ener       = araw.ener;
             temp       = araw.temp;
             runsfound0 = araw.n;
             
             if size(time,2) >= 201
             time0 = time(201:end) - time(201);
             ener0 = ener(2,201:end);
             
             x = time0';
             y = ener0';
             
             c = cESinf;
             a = cESini - cESinf;
             
             [fit1,fit2] = EfitconstraintES.fit(x,y,a,c);
             
             ESres{iangles,ibeta1,ibetaES,iVgs,1} = fit1;
             ESres{iangles,ibeta1,ibetaES,iVgs,2} = fit2;
             
             disp(['rsquare ', num2str(fit1.gof.rsquare), ...
                          ' ', num2str(fit2.gof.rsquare), ...
                    '  rmse ', num2str(fit1.gof.rmse), ...
                          ' ', num2str(fit2.gof.rmse)]);
             end
         end
      end
   end
end
    

save('data_ES/esFits_new.mat','nangles','beta1s','betaESs','Vgss','ESres');

%%

load('data_ES/esFits_new.mat');

%%
%a1 = cell(size(nangles,2),size(beta1s,2),size(betaESs,2),size(Vgss,2));
%a2 = cell(size(nangles,2),size(beta1s,2),size(betaESs,2),size(Vgss,2));
ndata = cell(size(nangles,2)*size(beta1s,2)*size(betaESs,2)*size(Vgss,2),4);
ic = 0;


for iangles = 1:size(nangles,2)
   for ibeta1 = 1:size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
             Clib.nangles = nangles(iangles);
             Clib.beta1   = beta1s(ibeta1);
             Clib.betaES  = betaESs(ibetaES);
             Clib.Vgs     = Vgss(iVgs);
             
             disp(['loading nangles ',num2str(Clib.nangles), ...
                 ' beta1 ', num2str(Clib.beta1), ...
                 ' betaES ', num2str(Clib.betaES), ...
                 ' Vgs ', num2str(Clib.Vgs)] );
             ic = ic + 1;
             a1 = ESres{iangles,ibeta1,ibetaES,iVgs,1};
             a2 = ESres{iangles,ibeta1,ibetaES,iVgs,2};
             
             ndata{ic,1} = iangles;
             ndata{ic,2} = ibeta1;
             ndata{ic,3} = ibetaES;
             ndata{ic,4} = iVgs;
                          
             b1(ic) = a1.gof.rsquare;
             b2(ic) = a2.gof.rsquare;
             
             disp(['rsquare ', num2str(a1.gof.rsquare), ...
                          ' ', num2str(a2.gof.rsquare), ...
                    '  rmse ', num2str(a1.gof.rmse), ...
                          ' ', num2str(a2.gof.rmse)]);
         end
      end
   end
end

%%

[bb1,ib1] = sort(b1);
[bb2,ib2] = sort(b2);


%% 
%New Plotting

for iangles = 9; %:13%size(nangles,2)
   for ibeta1 = 1; %:1%size(beta1s,2)
      for ibetaES = 3; %:1%size(betaESs,2)
         for iVgs = 1; %:3%size(Vgss,2)
             Clib.nangles = nangles(iangles);
             Clib.beta1   = beta1s(ibeta1);
             Clib.betaES  = betaESs(ibetaES);
             Clib.Vgs     = Vgss(iVgs);
             
             disp(['loading nangles ',num2str(Clib.nangles), ...
                 ' beta1 ', num2str(Clib.beta1), ...
                 ' betaES ', num2str(Clib.betaES), ...
                 ' Vgs ', num2str(Clib.Vgs)] );
             cESinf = esConst_inf{iangles,ibeta1,ibetaES,iVgs}.emean;
             cESini = esConst_ini{iangles,ibeta1,ibetaES,iVgs}.emean;
             araw   = rawES{iangles,ibeta1,ibetaES,iVgs};
             time       = araw.time;
             ener       = araw.ener;
             temp       = araw.temp;
             runsfound0 = araw.n;
             
             fit1 = ESres{iangles,ibeta1,ibetaES,iVgs,1};
             fit2 = ESres{iangles,ibeta1,ibetaES,iVgs,2};
             
             time0 = time(201:end) - time(201);
             ener0 = ener(2,201:end);
             
             x = time0';
             y = ener0';
             
             c = cESinf;
             a = cESini - cESinf;
             
             figure(100)
             plot(x,y,'b.');
             hold on
             %plot([x(1);x(end)],[c,c],'r');
             hold on
             %plot([x(1);x(end)],[a,a],'g');
             hold on
             %plot(fit1.fit);
             plot(fit2.fit,'b');           
         end
      end
   end
end

%%
plot(x,y,'ro')
xlabel('time (ps)');
ylabel('ES Energy (kcal/mol)' );
%%
plot(fit1.fit,'g.');
xlabel('time (ps)');
ylabel('ES Energy (kcal/mol)' );
%%
plot(fit2.fit,'b.');
xlabel('time (ps)');
ylabel('ES Energy (kcal/mol)' );
%% 
%Constants

ic = 0;
for iangles = 1:size(nangles,2)
   for ibeta1 = 1:size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
             Clib.nangles = nangles(iangles);
             Clib.beta1   = beta1s(ibeta1);
             Clib.betaES  = betaESs(ibetaES);
             Clib.Vgs     = Vgss(iVgs);
             
             ic = ic + 1;
             
             disp(['loading nangles ',num2str(Clib.nangles), ...
                 ' beta1 ', num2str(Clib.beta1), ...
                 ' betaES ', num2str(Clib.betaES), ...
                 ' Vgs ', num2str(Clib.Vgs)] );
             cESinf = esConstant{iangles,ibeta1,ibetaES,iVgs}.emean;
             cESini   = esConst_ini{iangles,ibeta1,ibetaES,iVgs}.emean;
             a_err   = esConst_ini{iangles,ibeta1,ibetaES,iVgs}.esd/100;
             araw    = rawES{iangles,ibeta1,ibetaES,iVgs};
             time       = araw.time;
             ener       = araw.ener;
             temp       = araw.temp;
             runsfound0 = araw.n;
             
             fit1 = ESres{iangles,ibeta1,ibetaES,iVgs,1};
             fit2 = ESres{iangles,ibeta1,ibetaES,iVgs,2};
             
             time0 = time(201:end) - time(201);
             ener0 = ener(2,201:end);
             
             x = time0';
             y = ener0';
             
             c = cESinf;
             a = cESini - cESinf;
             d(ic) = a - a_err;
             
             a1 = fit1.fit.a1;
             delta_a(ic) = abs(a1 - a);
                  
         end
      end
   end
end


%% 
% Plot ES energy vs. Chain Length

col = {'b','g','m','r','k','c'};
sym = {'o','x','^'};
ic = 0;
for iangles = 1:size(nangles,2)
   for ibeta1 = 1:1 %size(beta1s,2)
      for ibetaES = 1:size(betaESs,2)
         for iVgs = 1:size(Vgss,2)
             Clib.nangles = nangles(iangles);
             Clib.beta1   = beta1s(ibeta1);
             Clib.betaES  = betaESs(ibetaES);
             Clib.Vgs     = Vgss(iVgs);
             
             ic = ic + 1;
             
             disp(['loading nangles ',num2str(Clib.nangles), ...
                 ' beta1 ', num2str(Clib.beta1), ...
                 ' betaES ', num2str(Clib.betaES), ...
                 ' Vgs ', num2str(Clib.Vgs)] );
             cESinf = esConstant{iangles,ibeta1,ibetaES,iVgs}.emean;
             cESini   = esConst_ini{iangles,ibeta1,ibetaES,iVgs}.emean;
             araw    = rawES{iangles,ibeta1,ibetaES,iVgs};
             time       = araw.time;
             ener       = araw.ener;
             temp       = araw.temp;
             runsfound0 = araw.n;
             
             figure(301)
             hold on;
             plot(nangles(iangles),cESinf,[col{ibetaES},sym{iVgs}]);
         end
      end
   end
end

%% 
% tau of fit data

tau1 = zeros(size(nangles,2), size(betaESs,2), 3); % last is value, up, low
amp1 = zeros(size(nangles,2), size(betaESs,2), 3); % last is value, up, low
tau2 = zeros(size(nangles,2), size(betaESs,2), 3); % last is value, up, low
amp2 = zeros(size(nangles,2), size(betaESs,2), 3); % last is value, up, low

ibeta1 = 1; % beta1 = 1
iVgs   = 1; % Vgs = 0.6

for iangles = 1:size(nangles,2)
    for ibetaES = 1:size(betaESs,2)
       raw1 = rawES{iangles,ibeta1,ibetaES,iVgs};
       disp([num2str(iangles),' ',num2str(ibetaES)]);
       if (length(raw1.time) >= 201)
        fit1 = ESres{iangles,ibeta1,ibetaES,iVgs,1};
        fit2 = ESres{iangles,ibeta1,ibetaES,iVgs,2};
        ci  = confint(fit2.fit);
        
        t1 = max((ci(2,:)-ci(1,:))./(ci(2,:)+ci(1,:)));
        if ((t1 == Inf) | size( find(isnan(ci),1) ~= 0))
            t1 = 100;
        end
        maxerr= t1;
 
        if (maxerr < 0.5)
           disp(10);
            tau1(iangles, ibetaES, 1) = 1.0/fit2.fit.g21;
            tau1(iangles, ibetaES, 2) = abs(1.0/fit2.fit.g21 - 1.0/ci(2,2));
            tau1(iangles, ibetaES, 3) = abs(1.0/fit2.fit.g21 - 1.0/ci(1,2));
            tau2(iangles, ibetaES, 1) = 1.0/fit2.fit.g22;
            tau2(iangles, ibetaES, 2) = abs(1.0/fit2.fit.g22 - 1.0/ci(2,3));
            tau2(iangles, ibetaES, 3) = abs(1.0/fit2.fit.g22 - 1.0/ci(1,3));
%             amp1(iangles, ibetaES, 1) = fit2.fit.a21;
%             amp1(iangles, ibetaES, 2) = abs(fit2.fit.a21 - ci(1,1));
%             amp1(iangles, ibetaES, 3) = abs(fit2.fit.a21 - ci(2,1));
%             amp2(iangles, ibetaES, 1) = fit2.fit.a22;
%             amp2(iangles, ibetaES, 2) = 0.0;
%             amp2(iangles, ibetaES, 3) = 0.0;
            disp(11);
        else
           disp(12)
            ci = confint(fit1.fit);
            tau2(iangles, ibetaES, 1) = 1.0/fit1.fit.g1;
            tau2(iangles, ibetaES, 2) = abs(1.0/fit1.fit.g1 - 1.0/ci(2,1));
            tau2(iangles, ibetaES, 3) = abs(1.0/fit1.fit.g1 - 1.0/ci(1,1));
            tau1(iangles, ibetaES, 1) = 0.0;
            tau1(iangles, ibetaES, 2) = 0.0;
            tau1(iangles, ibetaES, 3) = 0.0;
%             amp2(iangles, ibetaES, 1) = fit1.fit.a1;
%             amp2(iangles, ibetaES, 2) = 0.0;
%             amp2(iangles, ibetaES, 3) = 0.0;
%             amp1(iangles, ibetaES, 1) = 0.0;
%             amp1(iangles, ibetaES, 2) = 0.0;
%             amp1(iangles, ibetaES, 3) = 0.0;
            disp(13)
        end
       end
    end
end

%%
% Plotting tau vs chain length

col = {'b','g','b','m','c','k','r','g','b','m','c','k','r'};
for i1 = 1:size(nangles,2)
    for i2 = 1:1%size(betaESs,2)
        figure(201)
        hold on;
        errorbar(nangles(i1),tau1(i1,i2,1),tau1(i1,i2,2),...
            tau1(i1,i2,3),[col{i2},'o']);
        figure(201)
        hold on;
        errorbar(nangles(i1),tau2(i1,i2,1),tau2(i1,i2,2),...
            tau2(i1,i2,3),[col{i2},'x']);
%         figure(301)
%         hold on;
%         errorbar(nangles(i1),amp1(i1,i2,1),amp1(i1,i2,2),...
%             amp1(i1,i2,3),[col{i2},'o']);
%         figure(302)
%         hold on;
%         errorbar(nangles(i1),amp2(i1,i2,1),amp2(i1,i2,2),...
%             amp2(i1,i2,3),[col{i2},'x']);
    end
end
figure(201)
title('tau 1');
% figure(202)
% title('tau 2');
% figure(301)
% title('amp 1');
% figure(302);
% title('amp 2');


