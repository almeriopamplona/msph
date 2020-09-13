%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  SMOOTHED PARTICLE HYDRODYNAMICS                    %%%
%%%               RUGE-KUTTA 2nd ORDER TIME INTEGRATION                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code integrates the properties of the particles through time using %
% the 2nd order Runge-Kutta method. The method is divided into an         %
% intermediate step and a final step. Therefore, besides the other inputs,%
% one should pay attention to "opt", which controls which step model to   %
% use.                                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% N        : Total particles number                        [int value   ] %            
% dt       : Time step                                     [double value] %     
% neighbor : Neighbors number and identification           [array       ] %
% part     : Properties of the particles                   [struct array] %    
% D        : Derivatives of the properties of the particles[struct array] %    
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% x        : Trajetory                                     [array]        %
% v        : Velocity                                      [array]        %
% e        : Internal energy                               [array]        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,e,x] = integration(opt,part,ipart,N,D,dt)  

  u = zeros(N,1);
  e = zeros(N,1);
  x = zeros(N,1);

if (opt == 1) % First step of the 2nd order Runge-Kutta
   for i = 1:N
     u(i) = ipart.u(i) + dt*D.u(i)/2; % Velocity   at t + dt/2
     e(i) = ipart.e(i) + dt*D.e(i)/2; % Energy     at t + dt/2
     x(i) = ipart.x(i) + dt*D.x(i)/2; % Position   at t + dt/2
   end
     
else          % Second step of the 2nd order Runge-Kutta
   for i = 1:N 
     u(i) = (ipart.u(i) + part.u(i) + dt*D.u(i))/2; % Velocity at t + dt
     e(i) = (ipart.e(i) + part.e(i) + dt*D.e(i))/2; % Energy   at t + dt 
     x(i) = (ipart.x(i) + part.x(i) + dt*D.x(i))/2; % Position at t + dt
   end
end
  
end
