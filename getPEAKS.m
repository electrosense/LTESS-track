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
  

% getPEAKS: get PSS peaks jumping from the last position known that was a valid peak
%
% Input:
%     - data: vector of IQ samples of th signal
%     - pss_step:  distance between PSS in I/Q samples
%     - search_window: range around peak to make analysis
%     - correlation_factor: value to determine if there is a peak or not
%     - reample_factor: resample factor to apply
%     - Z: Zadoof sequency in frequency
%     - Z_t: Zadoof sequency in time
%     - flip: (debug parameter) to flip the signal and zadoof sequencies
%     - SNR: value to add gaussian noise to the signal
%     
% Output:
%     - peaksV: vector with the correlation value of the peaks found
%     - peaksL: vector with the location value of the peaks found
%     - S: (debug) mean of the signal power (peaks)
%     - N: (debug) mean of the noise power.

function [peaksV, peaksL, S, N] = getPEAKS(data, pss_step, search_window, correlation_factor, resample_factor, Z, Z_t, flip, SNR)

    SAFE_WINDOW = 60;      
    
    index_done = 0;
    
    start_i = index_done + pss_step - search_window;    
    end_i = index_done + pss_step + search_window +1; 


    peak_list = [];
    
    S_values = [];
    N_values = [];
    
    mshort = Z_t;

    
    power_mshort=sqrt( sum(abs((mshort)).^2) );    

    myindex = 0;
    
    while (end_i< size(data,1))

        
        myindex = myindex + 1;
        slice = data(start_i:end_i);  
        % Add gaussian noise
        if (~ isnan (SNR))
           r_signal = awgn(real(slice), SNR,'measured');
           i_signal = awgn(imag(slice), SNR,'measured');
           slice = complex(r_signal, i_signal);
        end
        
        N =128;    % Size of mask
        
        xcor0 = (xcorr(slice,mshort));
        xcor0 = xcor0(size(slice):end);
        
        tmpnorm2 = sqrt(movsum((abs(slice)).^2, [1 N]));                    % power of local signal
        xcor_final = xcor0(1:size(tmpnorm2,1)) ./ (tmpnorm2*power_mshort);  % normalize the cross-corr
        
        res_corr = zeros(size(slice));
        res_corr(1:length(xcor_final)) = xcor_final;
                       
        % Upsampling if it is defined        
        if (resample_factor ~= 1)
            res_corr = resample(res_corr,resample_factor,1);
        end
                               
        res_corr = abs(res_corr);
                        
        % Find peaks  
        [p_value, p_loc] = findpeaks(res_corr,'MinPeakDistance',(size(slice,1)*resample_factor)-2,'MinPeakHeight', correlation_factor);  
  
        if (size(p_value,1) == 0)            
          fprintf('[WARNING] No peak detected, corr_factor=%.2f\n', correlation_factor);
        end
        
        if (size(p_loc,1) == 0)
            peak = NaN;
            p_value = NaN;
        else                   
            peak = start_i + (p_loc)/resample_factor -1;
        end        
        
        % Check if the peak is in the SAFE_WINDOW
        if( ~isnan(peak) && size(peak_list,1) >= 1)
            half = (start_i+end_i)/2.0;            
           
            if (abs ( half - peak )) > SAFE_WINDOW       
                %keyboard;
                fprintf('[WARNING] peak detected out of the SAFE_WINDOW(%d) -> %d\n',SAFE_WINDOW, abs(half-peak));                
                peak = NaN;
                p_value = NaN; 
            end
        end
        
        % Computing Signal and Noise terms (debug purposes ) 
        if (~ isnan(p_value))
           S_values = [S_values; p_value];
           N_values = [N_values; nanmean( res_corr( [1:p_loc-1 p_loc+1:end]))];
        end
        
        peak_list = [peak_list; p_value peak];
        
        if isnan(peak)
            % Current peak is nan -> No PSS detected
            % Move 9600 forward            
            start_i = start_i + pss_step;
            end_i = end_i + pss_step;            
        else           
            % Peak was found -> Move 9600 forward from here
            index_done = round(peak);
            start_i = index_done + pss_step - search_window;
            end_i = index_done + pss_step + search_window;    
        end
                
    end
    
    peaksV = peak_list(:,1);
    peaksL = peak_list(:,2);
    
    S = nanmean(S_values);
    N = nanmean(N_values);
    
    % Check if all the PEAKS are roughly 9600 IQ samples far between them.
    a = find(abs((peaksL(2:end)-peaksL(1:end-1))-pss_step) > 5 );
    
    if ( size(a,1) ~= 0)
        fprintf('[WARNING] Some peaks (%d) detected show a distance higher than 5 IQsamples between them\n',size(a,1));        
    end
    
end