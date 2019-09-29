function S = skew(v)
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% The function generates a skew-symmetric matrix for a give vector v
%	S = skew(v)

S = [0 -v(3) v(2); v(3) 0 -v(1);-v(2) v(1) 0];  
