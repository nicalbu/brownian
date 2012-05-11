classdef EfitES
   methods (Static)
      function y = testData(x, a21,t1, a22,t2, noise)
         [n1,n2] = size(x);
         y = a21 * exp(-x/t1) + a22 * exp(-x/t2) + noise*randn(n1,n2);
      end
         
      function [fit1struct, fit2struct] = fit(x,y,y0equil,yequil)         
         ft1 = fittype('a1*exp(-g1*x)+c',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a1', 'g1', 'c'});
         
         s1 = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[-Inf,0,-Inf],...
            'Upper',[Inf,Inf,Inf],...
            'Startpoint',[y0equil 1 yequil]);
         
         [fit1,gof1,out1] = fit(x,y,ft1,s1);
         fit1struct.x   = x;
         fit1struct.fit = fit1;
         fit1struct.gof = gof1;
         fit1struct.out = out1;
         
         ft2 = fittype('a2*((1-b)/2)*exp(-g21*x) + a2*((1+b)/2)*exp(-g22*x) +c',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a2', 'b', 'g21', 'g22', 'c'});
         
         s2 = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[-Inf,-1,0,0,-Inf],...
            'Upper',[Inf,1,Inf,Inf,Inf],...
            'Startpoint',[y0equil 0 1 1 yequil]);
         
         [fit2, gof2, out2] = fit(x,y,ft2,s2);
         fit2struct.x   = x;
         fit2struct.fit = fit2;
         fit2struct.gof = gof2;
         fit2struct.out = out2;
         
      end
      function y = ypred(fit1,x)
         y = fit1.fit(x);
      end
      function y = yobs(fit1)         
         y = Efit.ypred(fit1,fit1.x) - fit1.out.residuals;
         % tested using:
%          x = (1:100)';
%          y = Efit.testData(x, 1, 10, 1,20, 0.1);
%          [fit2, fit1] = Efit.fit(x,y);
%          ypred = Efit.ypred(fit2,x);
%          yobs  = Efit.yobs(fit2);
%          max(abs(y-yobs))
      end
      function [res, upper, lower] = value(fs, coefIn)
         % coef1 = a variable in the fit function (a1, g1, etc)
         %    if coef1 = t1, then it treats this an inverse of g1
         inverse = false;
         coef1 = coefIn;
         if (strcmp('t1',coefIn))
            inverse = true;
            coef1 = 'g21';
         end
         if (strcmp('t2',coefIn))
            inverse = true;
            coef1 = 'g22';
         end
         icoef = find(strcmp(coef1, coeffnames(fs.fit)));
         if (icoef == [])
            res = 0;
            lower = 0;
            upper = 0;
         else
            val = coeffvalues(fs.fit);
            ci = confint(fs.fit);
            if (~inverse)
               res = val(icoef);
               lower = res - ci(1,icoef);
               upper = ci(2,icoef) - res;
            else
               res = 1/val(icoef);
               lower = res - 1/ci(2,icoef);
               upper = 1/ci(1,icoef) - res;
            end
         end
      end
   end
end
