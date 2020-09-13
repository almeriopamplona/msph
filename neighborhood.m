%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  SMOOTHED PARTICLE HYDRODYNAMICS                    %%%
%%%                   NEAREST NEIGHBOUR PARTICLES                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     % 
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code finds the nearest neighbour particles of an ith particle,     %
% which belongs to a local or compact domain centred at the ith particle. % 
% The code uses a direct search method and builds a pseudo linked list    %
% with the total amount of particles in each compact domain and the index % 
% of each particle in these domains.                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %               %
% alpha    : SPH paramenters                            [array]           %
% N        : Total particles number                     [int value]       %                            
% part     : Particle position                          [struct array]    %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% neighbor : Neighbors number and identification        [array]           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neighbor] = neighborhood(part, N,alpha)

    neighbor = zeros(N,5); % Firts columm alocats the neighbors number and 
                           % the other columms alocat the neighbor ID.
    R     = N/alpha(4);    % Search radius
   
for i = 1:N
    for j = max(1,i-R):min(N,i+R) % Radius deslocation 
       if i~=j
          h = (part.h(i) + part.h(j))/2;
          if abs(part.x(i) - part.x(j)) <= alpha(5)*h % Suport domain
             neighbor(i,1) = neighbor(i,1) + 1; % Neighbors number
             neighbor(i,neighbor(i,1)+1) = j;   % Neighbors identification
          end
       end
    end
end
    
end
    
