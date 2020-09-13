%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              CORRECTED SMOOTHED PARTICLE HYDRODYNAMICS              %%%
%%%                          SMOOTHING KERNEL                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 27.12.2018                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% Subroutine to calculate :                                               %
%              * opt = 1 - cubic spline smoothing kernel                  %
%              * opt = 2 - gaussian smoothing kernel                      %
%              * opt = 3 - supergaussian smoothing kernel                 %
%              * opt = 4 - quartic smoothing kernel                       %
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

function f = W(opt,r,h)
  switch (opt)
% CUBIC SPLINE KERNEL --------------------------------------------------- % 

    case 1
      sigma = 1/h;
      s = r / h;
      if     s <= 1
        f = 2/3 - s^2 + 0.5*s^3;
      elseif s <= 2
        f = ((2 - s)^3)/6 ;
      else
        f = 0;
      end
      
% GAUSSIAN KERNEL ------------------------------------------------------- %

    case 2
      sigma = 1/ (h*sqrt(pi));
      s = r / h;
      if s <= 3
        f = exp(-s^2);
      else
        f = 0;
      end
      
% SUPERGAUSSIAN KERNEL -------------------------------------------------- %      
    
    case 3
      sigma = 1/(h*sqrt(pi));
      s = r / h;
      if s <= 3.0
        f = exp(-s^2) * (3/2 - s^2);
      else
        f = 0;
      end
      
% QUARTIC KERNEL -------------------------------------------------------- %      
    case 4
      sigma = 1/h;
      s     = r/h;
 
      if     (0 <= s)&&(s <= 2)
        f = (2/3) - (9/8)*s^2 + (19/24)*s^3 - (5/32)*s^4;
      elseif (s > 2)
        f = 0;
      end
  end
  f = f * sigma;
end
