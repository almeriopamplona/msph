%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                       HEAT CONDUCTION BALANCE                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code calculates the temperature derivative using the SPH approxi-  %
% mation of the energy balance.                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% h                : Smoothed lenght                             [double] %
% particle         : Particle position x-axis                    [struct] %
% numRealParticles : Real particles number                       [int]    %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% DT        : Lagrangian time derivative                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function DT = heat(opt,wsi,ws,numRealParticles,particle,neighbor,h,...
                   kernelDim,kernelOpt)

DT = zeros(numRealParticles,1);


for i = 1:numRealParticles
    for j = 1:neighbor(i,1)
        
      % neighbor counter        
        k    = neighbor(i,j+1); 
        
      % distance between particles on x-axis
        xji  = particle.x(k) - particle.x(i);
        
      % distance between particles on y-axis
        yji  = particle.y(k) - particle.y(i);
        
      % distance norm between particles 
        rji  = sqrt(xji^2 + yji^2);
                
      % Kernel gradient on x-axis
        dwxji = xji*DW(kernelDim,kernelOpt,rji,h)/(rji*h);            
        
      % Kernel gradient on y-axis
        dwyji = yji*DW(kernelDim,kernelOpt,rji,h)/(rji*h);
        
      % Kernel gradient vector
        dw  = [dwxji; dwyji];
        
        x = [xji; yji];
           
        if opt == 1
            
          % Renormalization matriz
            A   = [wsi(i,1) wsi(i,2); wsi(i,3) wsi(i,4)];      
            
          % Renormalizated kernel gradient vector
            dwc = A*dw;
            
        elseif opt == 2
            
          % Renormalization matriz
            B   = [-ws(i,2) -ws(i,3); -ws(i,4) -ws(i,5)]; 
            
          % Renormalizated kernel gradient vector
            dwc = B*dw;
            
         else
            
          % Nonrenormalizated kernel gradient
            dwc = dw;
            
         end
             
      % Corrected gradient kernel
        Fji = dot(x,dwc)/(rji^2 + (0.01*h)^2);
                
      % Conductivity component        
        if opt == 1
           kap = 4*particle.m*particle.ka^2/...
                (2*particle.rho^2*particle.Cp*particle.ka);
        elseif opt == 2
           kap = 4*particle.m*particle.ka^2/...
               (rji^2*particle.rho^2*particle.Cp*2*particle.ka);
        else
           kap = 4*particle.m*particle.ka^2/...
               (particle.rho^2*particle.Cp*2*particle.ka);
        end
        
      % heat conduction balance
        deltaT = particle.T(k) - particle.T(i);
               
        if opt == 1
            
            DT(i) = DT(i) + kap*( particle.T(i) - ...
                      particle.T(k) )*Fji*(1 - ws(i,6));
        else
            DT(i) = DT(i) - kap*deltaT*Fji;
        end
             
    end
end

end
