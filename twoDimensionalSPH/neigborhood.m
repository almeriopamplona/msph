%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                      NEIGHBORHOOD CALCULATION                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
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
%                                                                         %
% N               : Total particles number                       [int]    %
% h               : Particle smoothing length                    [double] %
% particle        : Particle position x-axis                     [struct] %                                           
% numRealParticles: number of real particles                     [int]    %                                                   
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% neighbor  : Neighborhood description                            [array] %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function neighbor = neigborhood(particle, h, numRealParticles, N)

% Initial neighbors 
  neighbor = zeros(numRealParticles, 21);
  
% Determination of the distance between the neighbor particles
  for i = 1:numRealParticles       % real particles loop
      for j = 1 : N % all  particles loop 
        if i ~= j
        
        % Norm of the distance between neighbor particles  
           rij = sqrt((particle.x(i) - particle.x(j))^2 + ...
                      (particle.y(i) - particle.y(j))^2);
        
          if rij <= 2*h
          
        % Neighbor number of the i particles
          neighbor(i, 1) = neighbor(i, 1) + 1;
          
          % Neighbor identification
          neighbor(i, neighbor(i, 1) + 1) = j;
        
          end
        end
      end
  end

end
