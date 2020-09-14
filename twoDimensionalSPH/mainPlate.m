%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                   HEAT FLOW ON A RETANGULAR PLATE                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc; tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEOMETRICAL PLATE DATA

  load('TE_correto_1.mat');

  Lx  = 0.1;               % plate length on x-axis                [m]
  Ly  = 0.1;               % plate length on y-axis                [m]
  esp = 6e-3/1.2;          % space between the particles           [m]
  nx  = round(Lx/esp);     % particles number on x-axis            [-]
  ny  = round(Ly/esp);     % particles number on y-axis            [-]
  h   = 1.2*esp;           % smoothed length                       [m]
  initialTemperature = ...
      [0.25 0.50 0.75 1];  % initial temperature 
                           %  - T1: botton                         [ºC]
                           %  - T2: right                          [ºC]
                           %  - T3: top                            [ºC]
                           %  - T4: left                           [ºC]
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL PROPERTIES

% Three different materials were selected in order to study the SPH method
% applied on thermal conduction problems, so one meust choose among them:
% 1 - carbon steel, 2 - Titanium, 3 - Vanadium.

  optMaterial =  1;
  
  beta = 0.1;                 % time step constant                 [-]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME SPECIFICATION

  t = 0.0;                    % Initial time                       [s] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
% PARTICLES CREATION

  [particle, N, numRealParticles, numGhostParticles] = ...
                                 creation(nx, ny, Lx, Ly, optMaterial);

  xi          = 9*h^(1/3);                % numerical height       [m]       
  particle.m  = particle.rho*Lx*Ly*xi/N;  % particle mass          [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY CONDITIONS

  particle = initial(numRealParticles, N, Lx, Ly, ...
                             particle, initialTemperature );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEIGHBORHOOD DETERMINATION

  neighbor = neigborhood(particle, h, numRealParticles, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME STEP    

  dt = beta*particle.rho*particle.Cp*h^2/particle.ka;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZATION
    
%  [ws,wsi] = renorm(2,r,part,neighbor,h,m,rho);       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KERNEL OPTION

  kernelOpt = 1; 
  kernelDim = 2;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME INTEGRATION 

  z = waitbar(0,'loading');

  err   = 10;   
  k     = 0;
  DT2   = zeros(numRealParticles,1);
  DT1   = zeros(numRealParticles,1);
  auxParticle = particle;

  while t < 115.89
     
      waitbar(t/115.89)
          
  % DISTRIBUTION PLOT 

    figure(2)
    plotPlate(particle, 1, 0, numRealParticles, initialTemperature, k, t);
    
  % IMPROVED EULER TIME INTEGRATION PART 1    
    oldParticle = particle;
    
  % NORMALIZATION 
    [ws,wsi] = renorm(2,numRealParticles,particle,neighbor,h, ...
                      kernelDim,kernelOpt); 
    
  % LAGRANGIAN TIME DERIVATE
    DT1 = heat(0,wsi,ws,numRealParticles,particle,neighbor,h, ...
               kernelDim,kernelOpt);     

  % EULERIAN TIME INTEGRATION
    
    %part(1:r,3) = integration(r,part,dt,DT);

    auxParticle.T(1:numRealParticles) = integration_r(1,particle,...
                                            numRealParticles,DT1,DT2,dt);
    
  % TIME EVOLUTION
    t = t + dt/2;
    
  % IMPROVED EULER TIME INTEGRATION PART 2
    
    [ws,wsi] = renorm(2,numRealParticles,particle,neighbor,h, ...
                      kernelDim,kernelOpt); 
    
  % LAGRANGIAN TIME DERIVATE
    
    DT2 = heat(0,wsi,ws,numRealParticles,auxParticle,neighbor,h, ...
               kernelDim,kernelOpt);   

  % EULERIAN TIME INTEGRATION
    particle.T(1:numRealParticles) = integration_r(2,particle,...
                                              numRealParticles,DT1,DT2,dt);

  % TIME EVOLUTION
    t = t + dt/2;
       
    err = norm(particle.T(1:numRealParticles) - ...
               oldParticle.T(1:numRealParticles));
     
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close(z);
tc = toc;
