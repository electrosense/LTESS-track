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

clear all
close all
addpath ('./common/')

SAMPLING_RATE = 1.92e6;
RESAMPLE_FACTOR = 60;
PSS_STEP = 9600;
SEARCH_WINDOW = 150;
CORRELATION_FACTOR = 0.1;
CORRELATION_FACTOR_PREAMBLE = 0.1;
PREAMBLE=20; % Number of analyzed PSS before start to jump
FLIP=0; % 0 disabled 1: enabled  (testing purposes)
POLYNOMIAL_DEGREE=1;

[Z, Z_t] = get_Zadoof();


filename='../PSS-Fabio/dataset/silver-armasuisse.bin';

capbuf = spec_load(filename);
chunk = capbuf(1:1.92e6*1);


% First pass
[PPM,PSS_percent,data,Y,p,PPM2] = getDrift(chunk,SAMPLING_RATE, Z, Z_t, PREAMBLE, ... 
    PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
sprintf('PPM: %f [%f] \nPSS detected: %f\n', PPM, PPM2, PSS_percent)

% Second pass
T_s = 1/1.92e6;
delta_f=(PPM*1e-6)*806e6;
Z_t_rotated = {};
Z_t_rotated{1} = Z_t{1}.*exp(-j*2*pi*T_s*delta_f*[1:length(Z_t{1})]');
Z_t_rotated{2} = Z_t{2}.*exp(-j*2*pi*T_s*delta_f*[1:length(Z_t{2})]');
Z_t_rotated{3} = Z_t{3}.*exp(-j*2*pi*T_s*delta_f*[1:length(Z_t{3})]');

[PPM,PSS_percent,data,Y,p,PPM2, th_learned, Z_sequence, last, p_loc20 ] = getDrift(chunk,SAMPLING_RATE, Z, Z_t_rotated, ... 
    PREAMBLE, PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
sprintf('PPM: %f [%f] \nPSS detected: %f\n', PPM, PPM2, PSS_percent)


t = PSS_STEP / SAMPLING_RATE;
IQ_2_ns = (1/(1.92e6))*1e9;

figure;
subplot(2,2,1)
plot(data(:,1)*t,data(:,2),'o'); hold on
plot(data(:,1)*t,Y,'-','LineWidth',3); 
title(strcat('Linear Regression ', ''));
ylabel('Cumm drift (IQ samples)');
xlabel('Time (seconds)');
legend({'Drift','Slope'},'Location','northwest');
grid on
mysetpic;

subplot(2,2,3)
plot(data(:,1)*t,data(:,2),'.'); hold on
plot(data(:,1)*t,Y,'-'); 
title('Linear Regression: Zoom in');
ylabel('Cumm drift (IQ samples)');
xlabel('Time (seconds)');
legend({'Drift','Slope'});

legend({'Drift','Slope'},'Location','northwest');
grid on
mysetpic;

%figure
subplot(2,2,2)
plot(data(2:end,1)*t,(diff(data(:,1))-1),'.-');     
title('PSS lost');
xlabel('Time (seconds)');
ylim([0 1])
mysetpic;

subplot(2,2,4)
plot(data(:,1)*t, (data(:,2)-Y)*IQ_2_ns,'.');
grid on
title('Residuals');
xlabel('Time (seconds)');
ylabel('Residuals (nanoseconds)');
mysetpic;


