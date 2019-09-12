function S = skew(w)
%
% The function generates a skew-symmetric matrix for a give vector v
%	S = skew(v)

S = [0 -v(3) v(2); v(3) 0 -v(1);-v(2) v(1) 0];  
