classdef SO3FunCBF < SO3Fun
% a class representing fibre symmetric ODFs
%
% Syntax
%   SO3F = SO3FunCBF(h,r,weights,psi)
%   SO3F = SO3FunCBF(fibre,'halfwidth',10*degree)
%
% Input
%  h       - @vector3d
%  r       - @vector3d
%  weights - double
%  psi     - @S2Kernel
%  fibre   - @fibre
%
% Output
%  SO3F - @SO3FunCBF
%
% Example
%   
%   f = fibre.beta(cs);
%   SO3F = SO3FunCBF(fibre,'halfwidth',10*degree)
%

  properties
    h
    r
    psi = S2DeLaValleePoussinKernel('halfwidth',10*degree);
    weights = 1;
    antipodal = false
  end

  properties (Dependent = true)
    SLeft
    SRight
    bandwidth % harmonic degree
  end
  
  methods
    function SO3F = SO3FunCBF(h,r,weights,psi,varargin)
      
      if nargin == 0, return;end
      
      SO3F.h = h;
      SO3F.r = r;
      SO3F.psi = psi;
      SO3F.weights = weights;
      SO3F.SS = getClass(varargin,'specimenSymmetry',specimenSymmetry);
                  
    end
    
    function SO3F = set.SRight(SO3F,SRight)
      if isa(SRight,'crystalSymmetry')
        SO3F.h = Miller(SO3F.h,SRight);
      else
        SO3F.h = vector3d(SO3F.h);
        SO3F.h.opt.sym = SRight;
      end
    end
    
    function SRight = get.SRight(SO3F)
      if isa(SO3F.h,'Miller')
        SRight = SO3F.h.CS;
      else
        try
          SRight = SO3F.h.opt.sym;
        catch
          SRight = specimenSymmetry;
        end
      end
    end
    
    function SO3F = set.SLeft(SO3F,SLeft)
      if isa(SLeft,'crystalSymmetry')
        SO3F.r = Miller(SO3F.r,SLeft);
      else
        SO3F.r = vector3d(SO3F.r);
        SO3F.r.opt.sym = SLeft;
      end      
    end
    
    function SLeft = get.SLeft(SO3F)
      if isa(SO3F.r,'Miller')
        SLeft = SO3F.r.CS;
      else
        try
          SLeft = SO3F.r.opt.sym;
        catch
          SLeft = specimenSymmetry;
        end
      end
    end
    
    function L = get.bandwidth(SO3F)
      L= SO3F.psi.bandwidth;
    end
    
    function SO3F = set.bandwidth(SO3F,L)
      SO3F.psi.bandwidth = L;
    end
        
  end 
  
  methods (Static = true)
  
    SO3F = example(varargin)
    
  end

end
