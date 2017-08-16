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


function [Z, Z_t] = get_Zadoof()
    % Initialize 
    Zseqfiles = {'./common/lte/25-Zadoff.it',
        './common/lte/29-Zadoff.it',
        './common/lte/34-Zadoff.it'};

    itload(Zseqfiles{1});
    Z{1} = seq; 
    itload(Zseqfiles{2});
    Z{2} = seq; 
    itload(Zseqfiles{3});
    Z{3} = seq;


    Z_t = {};
    
    for i=1:size(Z,2)
        seq = Z{i};        
        seqz = [seq(32:end); zeros(128-63,1); seq(1:31)];
        Z_t{i} = ifft(seqz(1:end-1), 128);        
    end
    
