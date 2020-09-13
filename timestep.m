%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                             TIME STEP                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code gets the time step for the time integral part of the SPH      %
% method. The method is based on the sound speed of the particles because %
% of the big role it plays in shock wave problems. The shock wave travels %
% through space-time with several different sounds of speed due to the    %
% continuously changing in pressure, density and temperature of the fluid.% 
% Therefore, the stability of the direct numerical time integration       %
% depends on the sound of the speed of the problem. However, for some     %
% cases, the time step remains constant along the process, meaning that   %
% this code may be unused.                                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% N    : Total particles number                                 [int]     %
% D    : derivatives os the properties of the particles         [struct]  %
% part : properties of the particles.                           [struct]  %
% alpha: SPH paramentet array                                   [array]   %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% dt   : Time step                                              [double]  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dt = timestep(part,N,alpha,D)

    for i = 1:N
        dtc = min( part.h(i)/(abs(part.c(i))+ 0.6*alpha(1)*abs(part.c(i))));
        dti = sqrt(min(part.h(i)/D.u(i)));
    end
    
    dt  = min(0.4*dtc,0.25*dti);
end
