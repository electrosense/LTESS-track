%   Copyright (C) Electrosense 2017
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/. 
% 
%   Authors: Roberto Calvo-Palomino <roberto.calvo [at] imdea [dot] org>
%            Fabio Ricciato <fabio.ricciato [at] fri.uni-lj [dot] si>
%


% analyzeDrift: Analyze the peaks of PSS to estimate and compute linear regression
%
% Input:
%    - peaks: vector of peaks position founded.
%    - pss_step: distance between PSS in I/Q samples
%
% Output:
%    - DATA(:,1) -> Number of the PSS (iteration)
%    - DATA(:,2) -> Cummulative drift
%    - Y: polyval output
% 
function [DATA, Y, p] = analyzeDrift (peaks, pss_step, degree)

    pss_detect = peaks(:);
    pss_detect = pss_detect(:) - pss_detect(1);     

    x = (1:1:size(pss_detect,1));
    t = [x(:) pss_detect];
    cum_drift = (t(:,2)) - ((t(:,1)-1)*pss_step);    

    data = [x(:) cum_drift];
    z= isnan(data(:,2));
    data ( z, :) = [];

    % Linear regresion
    p = polyfit(data(:,1), data(:,2), degree);
    y = polyval(p,data(:,1));

    DATA = data;
    Y = y;
       
end

