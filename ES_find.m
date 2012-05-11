clear classes;

Clib.type = 'ediffuse';
Clib.nangles = 2;
Clib.Vgs = 0.6;
Clib.betaES = -10;
Clib.beta1 = 1;

nangles = [2 3 4 5 6 7 8 9 10 15 20 30 50];
beta1s  = [1 30];
betaESs = [-10 -20 -30];
Vgss    = [0.6 1.5 2.5];

rawES = cell(size(nangles,2),size(beta1s,2),size(betaESs,2),size(Vgss,2),5);
 

%Nsave: 1 = angles, 2 = vels, 3= ener, 4 = wf, 5 = cent
Clib.time1 = 1; % time step will be beta1 * Clib.time1
Clib.nsave1  = [1 1 1 0 0];
Clib.nsteps1 = 200;

Clib.time2 = 0.1;
Clib.nsave2  = [1 1 1 0 0];
Clib.nsteps2 = 1501;

%% Retrieve data
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
                [time,ener,temp, runsFound0] = ediffuseGet(Clib);
                a.n = runsFound0;
                a.time = time;
                a.ener = ener;
                a.temp = temp;
                rawES{iangles,ibeta1,ibetaES,iVgs} = a;
            end
        end
    end
end

%% Save raw data
save('data_ES/rawES_final.mat','rawES');

%% Total nruns for each configuration, plotted as a histogram
load('data_ES/rawES_final.mat');
ic =0;
for iangles = 1:size(nangles,2)
    for ibeta1 = 1:size(beta1s,2)
        for ibetaES = 1:size(betaESs,2)
            for iVgs = 1:size(Vgss,2)
               ic = ic+1;
               ndata(ic) = rawES{iangles,ibeta1,ibetaES,iVgs}.n; 
            end
        end
    end
end
hist(ndata);

%%
% Some plots of the raw data

load('data_ES/rawES_final.mat');

ii = 0;
for iangles = 1:size(nangles,2)
    for ibeta1 = 2 %:size(beta1s,2)
        for ibetaES = 1 %:size(betaESs,2)
            for iVgs = 2 %:size(Vgss,2)
                disp(['loading nangles ',num2str(nangles(iangles)), ...
                    ' beta1 ', num2str(beta1s(ibeta1)), ...
                    ' betaES ', num2str(betaESs(ibetaES)), ...
                    ' Vgs ', num2str(Vgss(iVgs))] );
                ii = ii + 1;
                time0 = rawES{iangles,ibeta1,ibetaES,iVgs}.time;
                ener0 = rawES{iangles,ibeta1,ibetaES,iVgs}.ener;
                if size(time0,2) >= 201
                    x = time0(201:end) - time0(201);
                    y = ener0(2,201:end);
                    figure(ii)
                    plot(x,y,'b.')
                end
            end
        end
    end
end


