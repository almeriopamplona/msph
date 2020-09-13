%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  SMOOTHED PARTICLE HYDRODYNAMICS                    %%%
%%%                   MOMENTUM AND ENERGY EQUATIONS                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %                                          
% The code calculates the velocity, the acceleration and the internal     %
% energy variation using the SPH numerical approximation of the momentum  %
% and the energy balance.                                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% alpha   : Correction parameter                           [array]        %           
% neighbor: Neighbors number and identification            [array]        %
% N       : Total particles number                         [int value]    %  
% part    : Properties of the particles                    [struct array] %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% D       : Derivatives of the properties of the particles [struct array] %   
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = balancef(part,N,neighbor,alpha)
    
    D.u = zeros(N,1); % velocity  derivative
    D.e = zeros(N,1); % energy    derivative 
    D.x = zeros(N,1); % position  derivative

    for i = 1:N  
      for j = 1:neighbor(i,1)  
        k = neighbor(i,j+1);            % Neighbor term
        
        % Mean smoothed length between i and j particles
        h = 0.5*(part.h(i) + part.h(k));
        
        % Absolute space difference between i and j particles 
        x = abs(part.x(i) - part.x(k)); 
       
        % Space difference between i and j particles 
        xij = part.x(i) - part.x(k);
        
        % Velocity difference between i and j particle
        uij = part.u(i) - part.u(k);
        
        % Mean density between i and j particle
        rho = 0.5*(part.d(i) + part.d(k));
        
        % Mean sound speed between i and j particle
        cij = 0.5*(part.c(i) + part.c(k));
        
        % Smoothing kernel gradient
        dw = xij*DW(alpha(3),x,h)/(x*h);
        
        % Artificial viscosity
        
        if     uij*xij < 0 
            
            muij  = -(uij*xij/h)/((x/h)^2 + 0.1^2);
            
            IIij  = (alpha(1)*cij*muij + alpha(2)*muij^2)/rho;
            
        elseif uij*xij >= 0
            
            IIij  = 0;
            
        end         
        
        % Momentum balance
        D.u(i) = D.u(i) - part.m(k)*(part.p(i)/part.d(i)^2 + ...
                 part.p(k)/part.d(k)^2 + IIij)*dw;        
        
        % Energy balance
        D.e(i) = D.e(i) + 0.5*part.m(k)*(part.p(i)/part.d(i)^2 + ...
                 part.p(k)/part.d(k)^2 + IIij)*uij*dw;
 
        % Position derivative
        D.x(i) = part.u(i); 
             
      end
    end   
end
