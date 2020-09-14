%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                       PARTICLES INITIAL STATE                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code set the correct temperature on the boundaries of the          %
% rectangular plate. The side boundaries have continuous temperatures,    %
% while the edges have temperatures that follow a polar interpolation     %
% between the sides that form the edge.                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% N                  : Total particles number                    [int]    %
% Lx                 : Plate length on x-axis                    [double] %  
% Ly                 : Plate length on y-axis                    [double] %  
% numRealParticles   : Real particles number                     [int]    %
% initialTemperature : Initial temperatures                      [array]  %
% particle           : Properties of the particles               [struct] %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% particle           : Properties of the particlesv              [struct] %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function particle = initial( numRealParticles, N, Lx, Ly, ...
                             particle, initialTemperature )

% PLATE SIDE BORDERS TEMPERATURE 
  for i = (numRealParticles + 1):N    
    if (particle.y(i) > 0)
      if (particle.x(i) > Lx && particle.y(i) < Ly)
      % RIGHT SIDE ------------------------------------------------------ %
        particle.T(i) = initialTemperature(2);
      elseif (particle.y(i) > Ly && particle.x(i) > 0 && particle.x(i) < Lx)
      % TOP SIDE -------------------------------------------------------- % 
        particle.T(i) = initialTemperature(3);
      elseif (particle.x(i) < 0 && particle.y(i) < Ly)
      % LEFT SIDE ------------------------------------------------------- %
        particle.T(i) = initialTemperature(4);
      end
    else
    % BOTTOM SIDE ------------------------------------------------------- %
      if (particle.x(i) > 0 && particle.x(i) < Lx)
        particle.T(i) = initialTemperature(1);
      end
    end
  end

% PLATE EDGE TEMPERATURE  
 for i = numRealParticles + 1 : N
  % NORTHWEST ----------------------------------------------------------- %    
    if (particle.x(i) < 0 && particle.y(i) > Ly)
      theta = atan((particle.y(i) - Ly) / particle.x(i)) - pi;
      
      particle.T(i) = 2 * (initialTemperature(4) - initialTemperature(3)) /...
      pi * theta + (3 * initialTemperature(4) - 2 * initialTemperature(3));

  % NORTHEAST ----------------------------------------------------------- %
    elseif (particle.x(i) > Lx && particle.y(i) > Ly)
      theta = atan((particle.y(i) - Ly) / (particle.x(i) - Lx));

      particle.T(i) = 2 * (initialTemperature(3) - initialTemperature(2)) /...
      pi * theta + initialTemperature(2);
  
  % SOUTHWEST ----------------------------------------------------------- %
    elseif (particle.x(i) < 0 && particle.y(i) < 0)
      theta = atan(particle.y(i) / particle.x(i)) - pi;
   
      particle.T(i) = -2 * (initialTemperature(4) - initialTemperature(1)) /...
      pi * theta + (2 * initialTemperature(1) - initialTemperature(4));
   
  % SOUTHEAST ----------------------------------------------------------- % 
    elseif (particle.x(i) > Lx && particle.y(i) < 0)
      theta = atan(particle.y(i) / (particle.x(i) - Lx));
     
      particle.T(i) = -2 * (initialTemperature(1) - initialTemperature(2)) /...
      pi * theta + initialTemperature(2);
    end
 end
 
end
