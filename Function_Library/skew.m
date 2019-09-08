function R = skew(w)
% $Id: skew.m,v 1.1 2009-03-17 16:40:18 bradleyk Exp $
% Copyright (C) 2005, by Brad Kratochvil
%
%SKEW  generates a skew-symmetric matrix given a vector w
%
%	R = SKEW(w)
%
% See also: ROTAXIS, SKEWEXP, SKEWCOORDS.
%


  if 3 ~= size(w,1)
    error('SCREWS:skew','vector must be 3x1')
  end
  
  if isnumeric(w)
    R = zeros(3,3);
  end
  
  R(1,2) = -w(3);
  R(1,3) =  w(2);
  R(2,3) = -w(1);

  R(2,1) =  w(3);
  R(3,1) = -w(2);
  R(3,2) =  w(1);

%   R(1,1) = 0;
%   R(2,2) = 0;
%   R(3,3) = 0;
  
end
