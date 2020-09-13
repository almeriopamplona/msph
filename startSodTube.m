%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                            INITIAL VALUE                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code sets the properties of each particle that modes a             %
% one-dimensional Sod's tube. The tube, as a standard data, has 1.2 m of  %
% length with the zero of the modelling set in its centre.  The left side,%
% [-0.6, 0), has 320 particles, with initial velocity as 0.0 m/s, density %
% equals to 1.0 kg/m3, internal energy equals to 2.5 kJ/kg, pressure      %
% equals to 1.0 kPa, mass equals to 0.001875 kg and sound speed equal to  %
% 1.0 m/s. Meanwhile, the right side, [0, 0.6], has 80 particles with     %
% initial velocity as 0.0 m/s, density equal to 0.25 kg/m3, internal      %
% energy equals 1.795 kJ/kg, pressure equals to 0.1795 kPa, mass equals to%
% 0.001875 kg and sound speed equals to 0.845 m/s.                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
% N      : Total particles number                                         %
% dx1    : Left  spacial step                                             %
% dx2    : Right spacial step                                             %
% gamma  : Gas constant property                                          %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
% part.x : Particle position                                              %
% part.d : Particle density                                               %
% part.p : Particle pressure                                              %
% part.u : Particle velocity                                              %
% part.e : Particle internal energy                                       %
% part.s : Particle sound speed                                           %
% part.m : Particle mass                                                  %
% part.h : Particle smoothing length                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function part = preProcess(N,dx1,dx2,gamma)

    
for i = 1:N
   if i < ((N/5)*4 + 1)
      %part.x(i,1) = -(0.6 - dx1/2) + dx1*(i-1);    % Inomogeneos particle 
                                                    % position
      part.x(i,1) = -0.6 + dx1*(i-1);               % Homogeneos particle 
                                                    % position
      part.d(i,1) = 1.0;                              % Particle density
      part.e(i,1) = 2.5;                            % Particle Internal 
                                                    % Energy
      part.p(i,1) = (gamma - 1)*part.d(i)*part.e(i);% Particle pressure
      part.u(i,1) = 0.0;                            % Particle velocity
      part.c(i,1) = sqrt((gamma - 1)*part.e(i));    % Sound speed
      part.m(i,1) = dx1;                            % Particle mass
      part.h(i,1) = 8*dx1;                          % Smoothed length
    else
      %part.x(i,1) = dx2*(i-((N/5)*4 + 1) + 0.5);   % Inomogeneos particle
                                                    % position
      part.x(i,1) = 0.0 + dx2*(i - (N/5)*4);        % Homogeneos particle 
                                                    % positon
      part.d(i,1) = 0.25;                           % Particle density
      part.e(i,1) = 1.795;                          % Particle Internal 
                                                    % Energy
      part.p(i,1) = (gamma - 1)*part.d(i)*part.e(i);% Particle pressure
      part.u(i,1) = 0.0;                            % Particle velocity
      part.c(i,1) = sqrt((gamma - 1)*part.e(i));    % Sound speed
      part.m(i,1) = dx1;                            % Particle mass
      part.h(i,1) = 2*dx2;                          % Smoothed length
    end
end
end
