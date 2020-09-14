%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                        2D RENORMALIZATION                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code is an application of the kernel and gradient corrections pre- %
% sented by Bonet(1999), which intend to preserve momentum in the particle%
% formulation of SPH.                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% h               : Smoothed lenght                              [double] %               
% opt             : Choose between 2 normalization results       [int]    %                                 
% particle        : Particle position x-axis                     [struct] %            
% numRealParticles: Real particles number                        [int]    %
%                                                                         %
% OUTPUT ---------------------------------------------------------------- %
%                                                                         %
% ws  : Normalized smoothed kernel                                [array] %
% wsi : Normalized gradient of the smoothed kernel                [array] %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ws,wsi] = renorm(opt,numRealParticles,particle,neighbor,h,...
                           kernelDim,kernelOpt)

  ws  = zeros(numRealParticles,6);
  wsi = zeros(numRealParticles,4);

% Renormalization smooting kernel term calculation
  for i = 1:numRealParticles
    for j = 1:neighbor(i,1)
        
        k = neighbor(i,j+1);            % Neighbor term
      
      % distance between particles on x-axis
        xij  = particle.x(i) - particle.x(k); 
        
      % distance between particles on y-axis
        yij  = particle.y(i) - particle.y(k);
        
      % distance norm between particles 
        rij  = sqrt(xij^2 + yij^2);
        
      % Kernel gradient on x-axis
        dwxij = xij*DW(kernelDim,kernelOpt,rij,h)/(rij*h); 
        
      % Kernel gradient on y-axis
        dwyij = yij*DW(kernelDim,kernelOpt,rij,h)/(rij*h);
        
      % Smoothing kernel value 
        w = W(kernelDim,kernelOpt,rij,h);
        
      % Renormalization smoothing kernel term value
        ws(i,1) = ws(i,1) + particle.m*w/particle.rho;
        
        ws(i,6) = ws(i,6) - particle.m*(xij*dwxij + yij*dwyij)/particle.rho;
      
      % Renormalization smoothing kernel gradient term values
        
        if opt == 1
            
            ws(i,2) = ws(i,2) - particle.m*xij^2*dwxij/(2*rho);
        
            ws(i,3) = ws(i,3) - particle.m*yij^2*dwxij/(2*rho);
        
            ws(i,4) = ws(i,4) - particle.m*xij^2*dwyij/(2*rho);
        
            ws(i,5) = ws(i,5) - particle.m*yij^2*dwyij/(2*rho);
            
        elseif opt == 2
            
            ws(i,2) = ws(i,2) - particle.m*xij*dwxij/particle.rho;
        
            ws(i,3) = ws(i,3) - particle.m*yij*dwxij/particle.rho;
        
            ws(i,4) = ws(i,4) - particle.m*xij*dwyij/particle.rho;
        
            ws(i,5) = ws(i,5) - particle.m*yij*dwyij/particle.rho;
            
        end           
     end
    
    % Determinant of the gradient normalization matrix    
    D = ws(i,2)*ws(i,5) - ws(i,3)*ws(i,4);
    
    % Calculation of the elements of the inverse gradient norm matrix
    if D > 10^(-20)
        wsi(i,1) =  ws(i,5)/D;
        wsi(i,2) = -ws(i,3)/D;
        wsi(i,3) = -ws(i,4)/D;
        wsi(i,4) =  ws(i,2)/D;
    else
        wsi(i,1) =  1;
        wsi(i,2) =  0;
        wsi(i,3) =  0;
        wsi(i,4) =  1;
    end
    
  end
