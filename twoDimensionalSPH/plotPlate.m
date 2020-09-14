%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                   HEAT FLOW ON A RETANGULAR PLATE                   %%%
%%%                        PLOT OF THE RESULTS                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code plots the results of heat transport on a rectangular plate.   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% k                 : loop counter                              [int]     %
% numRealParticles  : total of the real particles               [int]     %
% initialTemperature: initial temperaute set in the plate       [array]   %
% particle          : Properties of the particles               [struct]  %
% plotOption        : Choose among the 3 types of plot ->                 %
%                      = 1 : 2D interpolated plot                         %
%                      = 2 : 3D interpolated plot                         %
%                      = 3 : 3D scattered plot                            %
% movieMaker        : Option of print each plot along the simulation such %
%                     that:                                               %
%                      = 1 : print                                        %
%                      = 0 : does not print                               %
%                                                                         %
% OUTPUT: --------------------------------------------------------------- %
%                                                                         %
% plot of the temperature in function of the position of the particles    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotPlate(particle, plotOption, movieMaker, numRealParticles, ...
                   initialTemperature, k, t)

    switch plotOption
        % 2D INTERPOLATED PLOT ------------------------------------------ %
        case 1
          % create a scatter plot and vary the circle's color
            scatter(particle.x(1:numRealParticles), ...
                    particle.y(1:numRealParticles),[],...
                    particle.T(1:numRealParticles),'filled'); 
          % interpolation of a scateter data:                                                         
            F = scatteredInterpolant(particle.x(1:numRealParticles),...
                                     particle.y(1:numRealParticles),...
                                     particle.T(1:numRealParticles));
            [X,Y] = meshgrid(particle.x(1:numRealParticles), ...
                             particle.y(1:numRealParticles));                  
          % T as a function of the intepolated scatter data
            interpTemperature = F(X,Y);
          % pseudocolor plot
            peseudoColorPlot = pcolor(X,Y,interpTemperature);
          % blurring the color's shape of a pcolor plot
            shading interp;
          % define the range of colors
            colormap jet
            colorbar
            peseudoColorPlot.EdgeColor = 'none';
            title(['Plate temperature distribution at t = ' ....
            sprintf('%0.4f',t) 's.'])
            xlabel('x (m)')
            ylabel('y (m)')
         
        % 3D INTERPOLATED PLOT ------------------------------------------ %    
        case 2
            
            leng = length(particle.x(1:numRealParticles));
            x = linspace(particle.x(1),particle.x(numRealParticles),leng);
            y = linspace(particle.y(1),particle.y(numRealParticles),leng);
            y = y';
            [X,Y] = meshgrid(x,y);
          % interpolation of a scateter data 
            F = scatteredInterpolant(partticle.x(1:numRealParticles),...
                                     particle.y(1:numRealParticles),...
                                     particle.T(1:numRealParticles));
          % T as a function of the intepolated scatter data
            interpTemperature = F(X,Y);
            meshc(X,Y,interpTemperature)
            colormap jet
            title(['Plate temperature distribution at t = ' ....
                    sprintf('%0.4f',t) 's.'])
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('Temperature (ºC)')
            axis([particle.x(1) particle.x(numRealParticles) ...
                  particle.y(1) particle.y(numRealParticles)...
                  0 max(initialTemperature)])
              
        % 3D SCATTERED PLOT --------------------------------------------- %      
        case 3    
            scatter3(particle.x(1:numRealParticles), ...
                     particle.y(1:numRealParticles), ...
                     partticle.T(1:numRealParticles),...
                     100*particle.T(1:numRealParticles)',...
                     particle.T(1:numRealParticles)','filled')
            %plot3(part(1:r,1), part(1:r,2), part(1:r,3), 'k*')
            title(['Plate temperature distribution at t = ' ....
            sprintf('%0.4f',t) 's.'])
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('Temperature (ºC)')
            axis([particle.x(1) particle.x(numRealParticles) ...
                  particle.y(1) particle.y(numRealParticles)...
                  0 max(initialTemperature)])
    end
   
    if movieMaker == 1
        
        k = k + 1;
    
        if k < 10
            print(['BSPH0000' sprintf('%d', k) '.png'], '-dpng');
        elseif (k >= 10)&&(k < 100)
            print(['BSPH000' sprintf('%d', k) '.png'], '-dpng');
        elseif (k >= 100)&&(k < 1000)
            print(['BSPH00' sprintf('%d', k) '.png'], '-dpng');
        else 
            print(['BSPH0' sprintf('%d', k) '.png'], '-dpng');
        end
        
    end
    
end
