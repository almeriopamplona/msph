%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  SMOOTHED PARTICLE HYDRODYNAMICS                    %%%
%%%                             DENSITY                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This method calculates the density of the particles based on the local  %
% amount of particles in the compact domain. This approach guarantees a   %
% mass balance.                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% N       : Total particles number                               [int]    %
% alpha   : Correction parameter                                 [array]  %
% neighbor: Neighbors number and identification                  [array]  %             
% part    : Properties of the particles                          [struct] %  
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% rho     : Smoothed density                                     [array]  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho = density(part,N,neighbor,alpha)
    
rho = zeros(N,1); % density vector pre-alocation

for i = 1:N  
  for j = 1:neighbor(i,1)  
      k = neighbor(i,j+1);  % Neighbor term
        
      % Mean smoothed length between i and j particles
      h = 0.5*(part.h(i) + part.h(k));
        
      % Absolute space difference between i and j particles 
      x = abs(part.x(i) - part.x(k)); 
               
      % Smoothing kernel
      w  = W(alpha(3),x,h);
      
      % Density sumation in SPH terms
      rho(i) = rho(i) + part.m(k)*w;  
  end
end
     rho = rho + 0.0855; % Density correction
end
