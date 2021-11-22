% analyze_20210823_air_v3  visualizes and analyzes Tian's 20210823 data 
% set for Kevin's gQD@SiO2 sample at air/coverslip interface.
%
% Haw Yang
% Princeton University
% August 24, 2021
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SET DATA-SPECIFIC VARIABLES    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%//////////////////////////////////////////////////////////////////////////
N_files = 5; % total number of files
filename_pattern = 'Interface_air%i.txt';
%//////////////////////////////////////////////////////////////////////////

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              LOAD DATA               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% pre-allocate space for data struct
  % filename
  % counts
  % parms
  % chisq
  % eval
data = struct('filename', cell(1, N_files), ...
              'counts', cell(1, N_files), ...
              'parms', cell(1,N_files), ...
              'chisq', cell(1,N_files), ...
              'eval', cell(1,N_files), ...
              'tau_avg', cell(1,N_files) );

%%% load instrument response data
IRF500 = load('IRF_laser.txt');
irfbase = IRF500(:,2);
%%% set the time base according to instrument response data
timebase = IRF500(:,1);
%%% load data and pack them into the struct array
for i=1:N_files
  data(i).filename = sprintf(filename_pattern,i);
  tmp = load(data(i).filename);
  data(i).counts = tmp(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              DATA FITTING            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = timebase;
instrument_response = irfbase;

%%%%% assign experiment-dependent parameters
%%% hand-assign valid data chanels, backeted by c1 and c2 TCSPC chanels

%//////////////////////////////////////////////////////////////////////////
c1 = 75; 
c1b = 375;
c2 = 3850;
%//////////////////////////////////////////////////////////////////////////

Nc = c2 - c1 + 1;           % number of channels containing valid photons
%%% instrument response function
irf = zeros(length(timebase),1);
irf(c1:c2) = instrument_response(c1:c2)/sum(instrument_response(c1:c2));
[~,ind] = max(irf);
t_zero = t(ind); % t-zero as indicated by instrument response

%%%%% assign initial gueses
    % The fitting quality is critically dependent on the initial guess, it
    % turns out. This algorithm, while employing simulated annealing to 
    % reduce the impacts of initial guesses to the fitting quality, is
    % still not perfect. Therefore, for different interfaces, e.g.,
    % solution or air, carefully constructed initial guesses are still
    % needed. For the initial guesses below, they are based on an average
    % of the amplitude and lifetime components from a small set of trial
    % data. Still, the final fitting results seem to depend strongly on the
    % amplitude initial guess. It may be because the sample size is too
    % small?

%//////////////////////////////////////////////////////////////////////////
t0 = 0.;             % t-zero shift
a0 = [0.28, 0.24];             % initial guess for amplitudes, air
tau0 = [0.17, 4.70, 35.76]; % initial guess for lifetimes, air
%//////////////////////////////////////////////////////////////////////////

%%% number of fitting variables
nv = length(tau0);  

%%%%%%%%%% Loop Through All Data Files %%%%%%%%%%

for i=1:N_files

  fprintf( 1, '\nProcessing %s\n', data(i).filename );
  
  %%% estimate background using c1 and c1b, chanel bounds for background
  bkg = mean( data(i).counts(c1:c1b) ); 
  ydata = data(i).counts(c1:c2);
  N = sum(ydata) - Nc*bkg; % number of signal photons
  
  %%% the objective function for fiminsearch(), chi-squared metric
  f = @(x,t,irf,bkg,N,c1,c2) ...
     mean( (model_trunc(t,irf,bkg,N,c1,c2, x(1),x(2:nv),x(nv+1:2*nv)) - ...
           ydata).^2./ydata );

  %%% function to be minimized, with additional parameters passed to it
  fun = @(x) f(x,t,irf,bkg,N,c1,c2);

  %%%%% do fitting by minimizing chi-squared metric 
  %%% set up initial guess for the fitting variables
  x0 = [t0, a0, tau0];

  %%% do simulated annealing to refine the initial guess
  fprintf( 1, ['  Performing simulated annealing for initial guess,' ...
               'simulannealbnd().\n'] );
  soptions = optimset('display','off');
  lb = [-1 0 0 0 1 5];
  ub = [1 0.8 0.8 1 5 50];
  [x,fval,exitFlag,output] = simulannealbnd(fun,x0,lb,ub,soptions);

  %%% do simplex minimization for final fitting
  fprintf( 1, ['  Performing simplex minimization for final fit,' ...
               'fminsearch().\n'] );
  options = optimset('MaxFunEvals',1e6,'MaxIter',1e6);
  x = fminsearch(fun,x,options);

  %%%%% pack fitting results to the data structure
  %%% the mean chi-squared metric as goodness of fit
  chisq = fun(x);
  %%% generate the model based on best-fit parameters
  y = model( t, irf, bkg, N, c1, c2, x(1), x(2:nv), x(nv+1:2*nv) );
  %%% set the time-delay axis based on best-fit parameter
  time_delay = t - t_zero - x(1); % time-delay t0 after fitting
  %%% sort the fitted (lifetime, weight) in increasing lifetime
  [~,ind] = sort(x(4:6)); % this only works for 3-component lifetime
  weight = [x(2:3) 1-sum(x(2:3))];
  lifetime = x(4:6);
  data(i).parms = [x(1) weight(ind) lifetime(ind)];
  data(i).chisq = chisq;
  data(i).eval = [time_delay; y];
  data(i).tau_avg = sum( weight.*lifetime );
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %         PLOT FITTING RESULTS         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  set(0, 'DefaultFigureRenderer', 'painters'); % vectorial EPS
  figure(1); clf;

  %%% overlay plotting data and fit
  subplot(4,1,1:3);
  semilogy( time_delay, instrument_response, 'Color', [0.5 0.5 0.5] ); 
  hold on;
  semilogy( time_delay, data(i).counts, 'b' );
  semilogy( time_delay, y, 'r', 'LineWidth', 2 ); 
  hold off;
  yy = ylim;
  axis([min(time_delay) max(time_delay) yy(1) yy(2)]);
  legend( 'IRF', 'data', 'fit', 'FontSize', 12 );
  set(gca,'XTick',[]);
  ylabel( 'photon counts' );
  set(gca, 'FontSize', 12 );
  title( 'analyze\_20210823\_air\_v3.m' );

  %%% plot the residuals
  subplot(4,1,4)
  plot( time_delay(c1:c2), ydata - y(c1:c2), 'k' );
  yy = ylim;
  axis([min(time_delay) max(time_delay) yy(1) yy(2)]);
  legend( ['\langle\chi^2\rangle' sprintf(' = %.2f',chisq)], ...
          'Location', 'NorthEast', 'FontSize', 12 );
  xlabel( 'time delay (ns)' );
  title( data(i).filename, 'Interpreter', 'none' );
  set(gca, 'FontSize', 12 );

  fprintf(1,'  All done.\n  Press any key to continue...\n' );
  drawnow;
  
  %%% save PNG figure to enable quick view of the fitting quality
  saveas(gcf,[data(i).filename '.png']);
  %%% save EPS figure to allow inclusion in a paper
  saveas(gcf,[data(i).filename '.eps'], 'epsc');

  % remove the commet "%" below if interactive visualization is desired.
  % pause; 
  
end

%%%%% preliminary analysis and statistics
lifetime = zeros(N_files,3);
weight = zeros(N_files,3);
tau_avg = zeros(N_files,1);
for i=1:N_files
  lifetime(i,:) = data(i).parms(5:7);
  weight(i,:) = data(i).parms(2:4);
  tau_avg(i) = data(i).tau_avg;
end
% %%% print out the sample average for weight and lifetime components
% mean(weight)
% mean(lifetime)
% %%% print out the sample mean of the lifetime
% mean( tau_avg )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         SAVE FITTING RESULTS         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% save as MATLAB binary
save analyze_20210823_air_v3.mat data;
%%% save as ASCII
fid = fopen( 'analyze_20210823_air_v3.out', 'w' );
fprintf(fid, '   t0 \t   a1 \t  tau1 \t   a2 \t  tau2 \t   a3 \t  tau3 \t chisq\n');
for i=1:N_files
  fprintf(fid, '%7.4f\t', data(i).parms(1)); % t0
  fprintf(fid, '%7.4f\t', data(i).parms(2)); % a1
  fprintf(fid, '%7.4f\t', data(i).parms(5)); % tau1
  fprintf(fid, '%7.4f\t', data(i).parms(3)); % a2
  fprintf(fid, '%7.4f\t', data(i).parms(6)); % tau2
  fprintf(fid, '%7.4f\t', data(i).parms(4)); % a3
  fprintf(fid, '%7.4f\t', data(i).parms(7)); % tau3
  fprintf(fid, '%5.2f\t', data(i).chisq);    % chisq
  fprintf(fid, '\n');
end
fclose(fid);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DEFINITION FOR IN-SCRIPT FUNCTIONS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = model( t, irf, bkg, N, c1, c2, t0, a, tau )
% MODEL generates lifetime model using the full extent of the time axis.
%  t, delay time array, parameter
%  irf, instrument response function array, parameter
%  bkg, background estimated from t<0, parameter
%  N, total number of signal photons
%  c1, beginning index for valid chanels
%  c2, ending index for valid chanels
%  t0, t-zero shift, variable
%  a, amplitude array, variable
%  tau, lifetime array, variable
  tt = t - t0;
  dt = t(2)-t(1);
  aa = [a, 1-sum(a)];
  y0 = 0;
  for i=1:length(tau)
    y0 = y0 + aa(i).*exp(-tt./tau(i))./tau(i).*dt;
  end
  y0 = y0 .* (tt >= 0);
  y1 = conv( y0, irf );
  y1 = y1 ./ sum(y1(c1:c2));
  y = y1(1:length(t)).*N + bkg;
end

function y = model_trunc( t, irf, bkg, N, c1, c2, t0, a, tau )
% MODEL_TRUNC generates lifetime model, truncated to valid-chanel range.
%  t, delay time array, parameter
%  irf, instrument response function array, parameter
%  bkg, background estimated from t<0, parameter
%  N, total number of signal photons
%  c1, beginning index for valid chanels
%  c2, ending index for valid chanels
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
  y1 = conv( y0, irf );
  y1 = y1 ./ sum(y1(c1:c2));
  y = y1(1:length(t)).*N + bkg;
  y = y(c1:c2);
end
