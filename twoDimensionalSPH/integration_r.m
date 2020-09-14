%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  SMOOTHED PARTICLE HYDRODYNAMICS                    %%%
%%%                    RUGE-KUTTA TIME INTEGRATION                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code calculates the time integration at t + dt/2 and at t + dt,    %
% using the second order Ruge-Kutta method.                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% dt       : Time step                                           [int]    %
% particle : Properties of the particles                         [struct] % 
% DT       : Lagrangian time derivative                          [array]  %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% T        : temperature distribution on the plate                [array] %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = integration_r(opt,particle,numRealParticles,DT1,DT2,dt)

if (opt == 1)
    
    T = particle.T(1:numRealParticles) + dt*DT1(1:numRealParticles);
    
else % Temperature at t + dt 
    
    T = particle.T(1:numRealParticles)+ 0.5*dt*(DT1(1:numRealParticles)...
        + DT2(1:numRealParticles));    
    
end
