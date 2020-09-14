%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                   PARTICLES GEOMETRY DETERMINATION                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code set the initial position of the particles, which forms a rec- %
% tangular plate.                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% nx        : Particles number on x-axis                         [int]    %  
% ny        : Particles number on y-axis                         [int]    %
% Lx        : Plate length on x-axis                             [double] %
% Ly        : Plate length on y-axis                             [double] %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% particle.x        : Particle position x-axis                   [struct] %              
% particle.y        : Particle position y-axis                   [struct] %            
% particle.T        : Particle temperature                       [struct] %          
% N                 : Total particles number                     [int]    %       
% numRealParticles  : Real particles number                      [int]    %      
% numGhostParticles : Ghost particles number                     [int]    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [particle, N, numRealParticles, numGhostParticles] = ...
                      rectangularPlateCreation(nx, ny, Lx, Ly,optMaterial)
 
  N     = (nx + 4) * (ny + 4);          % total particle number
  dx    = Lx/nx;
  dy    = Ly/ny;
  
  switch optMaterial
    
      case 1
          % Plain carbon steel
            particle.rho   = 7854;        % specific mass       [kg/m^3]
            particle.ka    = 60.5;        % thermal conductivit [W/(m*K)]
            particle.Cp    = 434;         % thermal capacity    [J/(kg*K)]
            particle.alpha = particle.ka/(particle.rho*particle.Cp); 
                                          % thermal difusivity  [m^2/s]

      case 2
          % Titanium
            particle.rho   = 4500;       % specific mass        [kg/m^3]
            particle.ka    = 21.9;       % thermal conductivit  [W/(m*K)]
            particle.Cp    = 522;        % thermal capacity     [J/(kg*K)]
            particle.alpha = particle.ka/(particle.rho*particle.Cp);
                                         % thermal difusivity   [m^2/s]

      case 3
          % Vanadium 
            particle.rho   = 6100;       % specific mass        [kg/m^3]
            particle.ka    = 30.7;       % thermal conductivit  [W/(m*K)]
            particle.Cp    = 489;        % thermal capacity     [J/(kg*K)]
            particle.alpha = particle.ka/(particle.rho*particle.Cp);
                                         % thermal difusivity   [m^2/s]

  end  
   
% Real particles  pre-allocation:  
  particleReal.x  = zeros(N,1);         % position x-axis          [m]
  particleReal.y  = zeros(N,1);         % position y-axis          [m]
  particleReal.T  = zeros(N,1);         % position temperature     [ºC]
% Ghost particles pre-allocation:  
  particleGhost.x = zeros(N,1);         % position x-axis          [m]
  particleGhost.y = zeros(N,1);         % position y-axis          [m]
  particleGhost.T = zeros(N,1);         % position temperature     [ºC]

% Number of particles:   
  numRealParticles  = 0;                % number of real particles 
  numGhostParticles = 0;                % number of ghost particles

% Properties allocation: 
  for i = 1:nx+4
    for j = 1:ny+4
      ix = dx*(i - 3 + 0.5);   % Index of the particles on x-axis
      jy = dy*(j - 3 + 0.5);   % Index of the particles on y-axis
      
      if (ix > 0 && ix < Lx && jy > 0 && jy < Ly) 
      % real particles
        numRealParticles = numRealParticles + 1;
        particleReal.x(numRealParticles) = dx*(i - 3 + 0.5); 
        particleReal.y(numRealParticles) = dy*(j - 3 + 0.5); 
        particleReal.T(numRealParticles) = 0;                
      else
      % phanton particles
        numGhostParticles = numGhostParticles + 1;
        particleGhost.x(numGhostParticles) = dx*(i - 3 + 0.5); 
        particleGhost.y(numGhostParticles) = dy*(j - 3 + 0.5); 
        particleGhost.T(numGhostParticles) = 0;                
      end
      
    end
  end
  
  particle.x = [ particleReal.x( 1:numRealParticles ) ; ...
                 particleGhost.x(1:numGhostParticles)]; 
  particle.y = [ particleReal.y( 1:numRealParticles ) ; ...
                 particleGhost.y(1:numGhostParticles)]; 
  particle.T = [ particleReal.T( 1:numRealParticles ) ; ...
                 particleGhost.T(1:numGhostParticles)];            
end
