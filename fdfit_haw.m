
%cw is the center postion of the instrument respons function
%y is the fluoresence decay of the sample
%how many exponentials you want in this fitting, 1==1exp 2==2exp 3==3exp.
%setpsize is the time interval between pixels.
close all
a=61;
b=3786;
t=aa_timeaxis(a:b,1);


yraw=aa_interface_airfix(a:b,1);
irfraw=aa_irf(a:b,2);

ydata=yraw./max(yraw);
irf=irfraw./sum(irfraw);
[~,ind] = max(irf);
t_zero = t(ind);
%//////////////////////////////////////////////////////////////////////////
t0 = 1;             % t-zero shift


   
           
a0 = [0.6, 0.2, 0.1];            
tau0 = [1, 19,32, 300];
bkg = 0.0001;
% initial guess for lifetimes, air
%//////////////////////////////////////////////////////////////////////////

x0 = [t0, a0, tau0, bkg];
 %%% the objective function for fiminsearch(), chi-squared metric

   f = @(x,t,irf) ...
     mean( (model_trunc(t,irf,x(9), x(1),x(2:4),x(5:8)) - ...
           ydata).^2./ydata );

  %%% function to be minimized, with additional parameters passed to it
  fun = @(x) f(x,t,irf);

  %%%%% do fitting by minimizing chi-squared metric 
  %%% set up initial guess for the fitting variables
  
  
 
lb = [0 0 0 1  1 10 -1];
   ub = [ 1 1 1 10 10 500 1000 1];

   
   y=fun(x0);
    soptions = optimset('display','off');
   [x,fval,exitFlag,output] = simulannealbnd(fun,x0,lb,ub,soptions);

  %%% do simplex minimization for final fitting
  fprintf( 1, ['  Performing simplex minimization for final fit,' ...
               'fminsearch().\n'] );
  options = optimset('MaxFunEvals',1e6,'MaxIter',1e6);
  x = fminsearch(fun,x,options);
  
  chisq = fun(x)*max(yraw);

 weight(1,:) = [x(2:4) 1-sum(x(2:4))];
 lifetime(1,:) = x(5:8);
 
 
 



plot(model_trunc(t,irf,x(9), x(1),x(2:4),x(5:8)))
hold on
plot(ydata)
set(gca, 'YScale', 'log');







function y = model_trunc( t, irf, bkg, t0, a, tau )
% MODEL_TRUNC generates lifetime model, truncated to valid-chanel range.
%  t, delay time array, parameter
%  irf, instrument response function array, parameter
%  bkg, background estimated from t<0, parameter
%  t0, t-zero shift, variable
%  a, amplitude array, variable
%  tau, lifetime array, variable
  tt = t - t0;
  dt = t(2)-t(1);
  y0 = 0;
  aa = [a, 1-sum(a)];
  for i=1:length(tau)
    y0 = y0 + aa(i).*exp(-tt./tau(i))./tau(i).*dt;
  end
  y0 = y0 .* (tt >= 0);
  y1 = conv( y0, irf )+bkg;
  y1 = y1 ./ max(y1);
  y = y1(1:length(t));
end


