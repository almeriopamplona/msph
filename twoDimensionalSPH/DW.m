%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   SMOOTHED PARTICLE HYDRODYNAMICS                   %%%
%%%                     SMOOTHING KERNEL GRADIENT                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Almério José Venâncio Pains Soares Pamplona                     %
% Date  : 29.06.2019                                                      %
% E-mail: almeriopamplona@gmail.com                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:                                                            %
%                                                                         %
% This code calculates the local smoothed gradient for any dimension (d = %
% 1, 2 or 3). The options of smoothed kernels are:                           %
%   - 1: cubic spline (Monaghan, 1992)                                    %
%   - 2: Gaussian (Liu, 2010)                                             %
%   - 3: Super Gaussian (Monaghan, 1992)                                  %
%   - 4: quintic spline (Liu, 2010)                                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                                  %
%                                                                         %
% r : Absolute value of the distance between two particles        [double]%    
% h : Average value of the smoothing length between two particles [double]%
%                                                                         %
% OUTPUT ---------------------------------------------------------------- %
%                                                                         %
% f : Smoothing kernel gradient result                            [double]%    
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = DW(d,opt,r,h)

% 1D SMOOTHED GRADIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if d == 1
  switch (opt)
      
  % CUBIC SPLINE GRADIENT (MONAGHAN, 1992) ------------------------------ %    
    case 1
      sigma = 2.0/(3.0*h);
      s = r / h;
      if ( s >= 0.0 &&  s <= 1.0 )
        f = -3.0*s + 9.0*0.25*s*s;
      elseif ( s > 1.0 && s <= 2.0 )
        f = -3.0*0.25*( (2.0 - s)^2 ) ;
      else
        f = 0.0;
      end
   
   % GAUSSIAN GRADIENT (LIU, 2010) -------------------------------------- %   
     case 2
        sigma = 1.0/(h*pi^0.5);
        s = r / h;
        if( s >= 0.0 && s <= 3.0  )
            f = -2.0*s*exp(-s*s);
        else
            f = 0.0;
        end
            
   % SUPERGAUSSIAN GRADIENT (MONAGHAN, 1992) ---------------------------- %       
     case 3
        sigma = 1.0/(pi^(0.5*d)*h^d);
        s = r / h;
        if( s >= 0.0 && s <= 3.0 )
            f = -2.0*exp(-s*s)*(d + 4 - 2.0*s*s);
        else
            f = 0.0;
        end
        
   % QUINTIC GRADIENT (LIU, 2010) --------------------------------------- %      
     case 4    
        sigma = 1.0/(120.0*h);
        s = r / h;
        if( s >= 0.0 && s <= 1.0)
            f = -50*s*s*s*s + 120*s*s*s - 120*s;
        elseif( s > 1.0 && s <= 2.0 )
            f = -5*(3.0 - s)^5 + 30*(2.0 - s)^4;
        elseif( s > 2.0 && s <= 3.0)
            f = -5*(3.0 - s)^4;
        else
            f = 0.0;
        end
  end
  
% 2D SMOOTHED GRADIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif d == 2
  switch (opt)
      
  % CUBIC SPLINE GRADIENT (MONAGHAN, 1992) ------------------------------ %    
    case 1
      sigma = 10.0/(7.0*pi*h^2);
      s = r / h;
      if ( s >= 0.0 &&  s <= 1.0 )
        f = -3.0*s + 9.0*0.25*s*s;
      elseif ( s > 1.0 && s <= 2.0 )
        f = -3.0*0.25*( (2.0 - s)^2 ) ;
      else
        f = 0.0;
      end
   
   % GAUSSIAN GRADIENT (LIU, 2010) -------------------------------------- %   
     case 2
        sigma = 1.0/(pi*h*h);
        s = r / h;
        if( s >= 0.0 && s <= 3.0  )
            f = -2.0*s*exp(-s*s);
        else
            f = 0.0;
        end
            
   % SUPERGAUSSIAN GRADIENT (MONAGHAN, 1992) ---------------------------- %       
     case 3
        sigma = 1.0/(pi^(0.5*d)*h^d);
        s = r / h;
        if( s >= 0.0 && s <= 3.0 )
            f = -2.0*exp(-s*s)*(d + 4 - 2.0*s*s);
        else
            f = 0.0;
        end
        
   % QUINTIC GRADIENT (LIU, 2010) --------------------------------------- %      
     case 4    
        sigma = 1.0/(478.0*pi*h*h);
        s = r / h;
        if( s >= 0.0 && s <= 1.0)
            f = -50*s*s*s*s + 120*s*s*s - 120*s;
        elseif( s > 1.0 && s <= 2.0 )
            f = -5*(3.0 - s)^5 + 30*(2.0 - s)^4;
        elseif( s > 2.0 && s <= 3.0)
            f = -5*(3.0 - s)^4;
        else
            f = 0.0;
        end
  end    

% 3D SMOOTHED GADIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif d == 3
  switch (opt)
      
  % CUBIC SPLINE GRADIENT (MONAGHAN, 1992) ------------------------------ %    
    case 1
      sigma = 1.0/(pi*h*h*h);
      s = r / h;
      if ( s >= 0.0 &&  s <= 1.0 )
        f = -3.0*s + 9.0*0.25*s*s;
      elseif ( s > 1.0 && s <= 2.0 )
        f = -3.0*0.25*( (2.0 - s)^2 ) ;
      else
        f = 0.0;
      end
   
   % GAUSSIAN GRADIENT (LIU, 2010) -------------------------------------- %   
     case 2
        sigma = 1.0/(h*h*h*pi^1.5);
        s = r / h;
        if( s >= 0.0 && s <= 3.0  )
            f = -2.0*s*exp(-s*s);
        else
            f = 0.0;
        end
            
   % SUPERGAUSSIAN GRADIENT (MONAGHAN, 1992) ---------------------------- %       
     case 3
        sigma = 1.0/(pi^(0.5*d)*h^d);
        s = r / h;
        if( s >= 0.0 && s <= 3.0 )
            f = -2.0*exp(-s*s)*(d + 4 - 2.0*s*s);
        else
            f = 0.0;
        end
        
   % QUINTIC GRADIENT (LIU, 2010) --------------------------------------- %      
     case 4    
        sigma = 1.0/(120.0*pi*h*h*h);
        s = r / h;
        if( s >= 0.0 && s <= 1.0)
            f = -50*s*s*s*s + 120*s*s*s - 120*s;
        elseif( s > 1.0 && s <= 2.0 )
            f = -5*(3.0 - s)^5 + 30*(2.0 - s)^4;
        elseif( s > 2.0 && s <= 3.0)
            f = -5*(3.0 - s)^4;
        else
            f = 0.0;
        end
  end
end 

  f = f*sigma;
end
