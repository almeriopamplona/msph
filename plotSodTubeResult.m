%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SMOOTHED PARTICLE HYDRODYNAMICS                      %
%                        SOD'S TUBE - RESULT PLOT                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     % 
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code is pos-processing plot of the numerical analysis with SPH ap- %
% plied on the Sod's tube case.                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSodTubeResults(part, data, makeVideo)

figure(1)
   
subplot(2,2,1)  % density plot 
plot(part.x,part.d,'ok',data.x,data.rho,'k-','LineWidth',...
     0.75,'MarkerSize',1.2)
%title(['t=' sprintf('%.4f', t) ' s.'])
xlabel('Position (m)')
ylabel('Density (kg/m^3)')
axis([-0.6 0.6 0.15 1.1])
   
subplot(2,2,2)  % pressure plot
plot(part.x,part.p,'ok',data.x,data.P,'k-','LineWidth',...
     0.75,'MarkerSize',1.2)
%title(['t=' sprintf('%.4f', t) ' s.'])
xlabel('Position (m)')
ylabel('Pressure (kPa)')
axis([-0.6 0.6 0.1 1.1])
  
subplot(2,2,3)  % velocity plot
plot(part.x,part.u,'ok',data.x,data.u,'k-','LineWidth',...
     0.75,'MarkerSize',1.2)
%title(['t=' sprintf('%.4f', t) ' s.'])
xlabel('Position (m)')
ylabel('Velocity (m/s)')
axis([-0.6 0.6 -0.05 0.9])
   
subplot(2,2,4)  % energy plot
plot(part.x,part.e,'ok',data.x,data.e,'k-','LineWidth',...
     0.75,'MarkerSize',1.2)
%title(['t=' sprintf('%.4f', t) ' s.'])
xlabel('Position (m)')
ylabel('Internal Energy (kJ/kg)')
axis([-0.6 0.6 1.7 2.7])

% Sequence print plot for video make ------------------------------------- %
if makeVideo == 1
    if k < 10
        print(['ASPH0000' sprintf('%d', k) '.png'], '-dpng');
    elseif (k >= 10)&&(k < 100)
        print(['ASPH000' sprintf('%d', k) '.png'], '-dpng');
    elseif (k >= 100)&&(k < 1000)
        print(['ASPH00' sprintf('%d', k) '.png'], '-dpng');
    else 
        print(['ASPH0' sprintf('%d', k) '.png'], '-dpng');
    end
end
% ------------------------------------------------------------------------ %
end
