%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                     SHOCK TUBE CONTROL PROGRAM                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona 
% Date  : 27.12.2018
% E-mail: almeriopamplona@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL DATA

n        = 80;                       % Minimum particle number
N        = 5*n;                      % Total particle number
dx1      = 0.6/(4*n);                % Left  spacial step  
dx2      = 0.6/n;                    % Right spacial step
gamma    = 1.4;                      % Gas constant property         
part     = start(N,dx1,dx2,gamma);   % Particle initial value
T        = 0.24;                     % Final   time
t        = 0.0;                      % Initial time
k        = 0.0;                      % Initial interaction count   
load('data_2.mat')                   % Analitical Solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALPHA PARAMETER VECTOR

% alpha(1) : alpha artificial viscosity parameter
  alpha(1)  = 1;
% alpha(2) : beta  artificial viscosity parameter
  alpha(2)  = 2;
% alpha(3) : kernel option: 
%               * 1 - cubic spline
%               * 2 - gaussian 
%               * 3 - supergaussinan 
%               * 4 - new quartic spline
  alpha(3)  = 4;
% alpha(4) : radius length of the neighborhood function (N/alpha(4)),
%            please, select a multiple of 5
  alpha(4)  = 20;
% alpha(5) : multiple of the smoothed length (h = alpha(6)*h)
%               * 2 - cubic spline
%               * 2 - gaussian
%               * 3 - supergaussian
%               * 2 - new quartic spline
  alpha(5)  = 2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL SOLUTION - SPH SCHEME WITH RUGE-KUTTA TIME INTEGRATION

z    = waitbar(0,'Loading'); % Waiting bar
D.u  = ones(N,1);

while t < T
    
    waitbar(t/T)
       
%%%%%%% TIME INTEGRATION FIRST STEP (RK2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Particle i neighborhood determination 
        neighbor = neighborhood(part,N,alpha);
    
        % Time step determination
        dt = timestep(part,N,alpha,D);
        
        % Density determination
        rho = density(part,N,neighbor,alpha);
        
        % Data actualization at t + dt/2
        ipart   = part;
        ipart.d = rho;
        ipart.p = (gamma-1)*ipart.d.*part.e;  % Pressure          [kPa]
        ipart.c = sqrt((gamma-1)*part.e);     % Sound speed       [m/s]
        
        % Trajetory, momentum and energy equations
        D = balancef(part,N,neighbor,alpha);
        
        % Time integration
        [v,e,x] = integration(1,part,ipart,N,D,dt);
        
        % Data actualization at t + dt/2
        ipart.x = x;                           % Position          [m]    
        ipart.u = v;                           % Velocity          [m/s]
        ipart.e = e;                           % Internal energy   [kJ/kg]
        
%%%%%%% TIME INTEGRATION SECOND STEP (RK2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        t = t + dt/2;   % first time actualization 
        
        % Particle i neighborhood determination 
        [neighbor] = neighborhood(ipart,N, alpha);
        
        
        % Density determination
        rho = density(ipart,N,neighbor,alpha);
        
        % Data actualization at t + dt
        part.d = rho;
        part.p = (gamma-1)*part.d.*ipart.e;  % Pressure          [kPa]
        part.c = sqrt((gamma-1)*ipart.e);    % Sound speed       [m/s]
        
        
        % Trajetory, continuity, momentum and energy equations
        D = balancef(ipart,N,neighbor,alpha);
        
        % Time integration
        [v,e,x] = integration(2,part,ipart,N,D,dt);
                        
        % Data actualization at t + dt
        part.x = x;                           % Position          [m]    
        part.u = v;                           % Velocity          [m/s]
        part.e = e;                           % Internal energy   [kJ/kg]
        
%%%%%%% Boundary condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q = 170;                  % number of particles at the left part, 
                                  % count from the first one particle
        m = 30;                   % number of particles at the right part, 
                                  % count from the last one particle
        part.d(1:q)      = 1;     % left  density  actualization
        part.p(1:q)      = 1;     % left  pressure actualization
        part.u(1:q)      = 0;     % left  velocity actualization
        part.e(1:q)      = 2.5;   % left  energy   actualization
        part.d((N-m):N)  = 0.25;  % right density  actualization
        part.p((N-m):N)  = 0.1795;% right pressure actualization
        part.u((N-m):N)  = 0;     % rigth velocity actualization
        part.e((N-m):N)  = 1.795; % rigth energy   actualization
         
%%%%%%% Error detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if (sum(sum(isnan(part.d))) > 0)||(sum(sum(isnan(part.p))) > 0)||...
   (sum(sum(isnan(part.u))) > 0)||(sum(sum(isnan(part.e))) > 0)||...
   (sum(sum(isnan(part.c))) > 0)
        
error('Error! NaN value detected (Interation %d, time: %.3f)\n\n', k, t);

        
else
            
   t = t + dt / 2; % second time actualization 
   k = k + 1;      % loop counter actualization 
        
end     


%%%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close(z); toc; % numerical time display

plotSodTubeResults(part, data, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
close(z); toc; % numerical time display
% FINAL DATA PRINT/SAVE

% save('teste1','part')                   % data save 
