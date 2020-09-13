%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              CORRECTED SMOOTHED PARTICLE HYDRODYNAMICS              %%%
%%%                     SMOOTHING KERNEL GRADIENT                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            % 
%                                                                         %
% Subroutine to calculate :                                               %
%              * opt = 1 - cubic spline smoothing kernel gradient         %
%              * opt = 2 - gaussian smoothing kernel gradient             %
%              * opt = 3 - supergaussian smoothing kernel gradient        %
%              * opt = 4 - quartic smoothing kernel gradient              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% r : Absolute value of the distance between two particles                %
% h : Average value of the smoothing length between two particles         %
%                                                                         %
% OUTPUT ---------------------------------------------------------------- %
%                                                                         %
% f : Smoothing kernel gradient result                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = DW(opt,r,h)
  switch opt
  
% CUBIC SPLINE KERNEL GRADIENT ------------------------------------------ %

    case 1
      sigma = 1/h;
      s = r/h;
      if     s <= 1
        f = - 2*s + 1.5*s^2;
      elseif s <= 2.0
        f = - 0.5*(2.0 - s)^2;
      else
        f = 0;
      end
       
% GAUSSIAN KERNEL GRADIENT ---------------------------------------------- %

    case 2
      sigma = 1/(h*sqrt(pi));
      s = r/h;
      if s <= 3
        f = - 2*s*exp(-s^2);
      else
        f = 0;
      end

% SUERGAUSSIAN KERNEL GRADIENT ------------------------------------------ %

    case 3
      sigma = 1/(h*sqrt(pi));
      s = r/h;
      if s <= 3
        f = (-2*s*(3/2 - s^2)*exp(-s^2)) - 2*s*exp(-s^2);
      else
        f = 0;
      end

% QUARTIC KERNEL GRADIENT ----------------------------------------------- %

    case 4
      sigma = 1/h;
      s     = r/h;
      if     (0 <= s)&&(s <= 2)
        f = - (9/4)*s + (19/8)*s^2 - (5/8)*s^3;
      elseif (s > 2)
        f = 0;
      end
  end
  f = f*sigma;
end
