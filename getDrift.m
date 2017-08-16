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

% getDrift: Compute the drift based on the PSS signal
% 
% Input:
%     - capbuf: vector of IQ samples
%     - sampling_rate: in Hz
%     - Z: Zadoof sequencies in frequency
%     - Z_t: Zadoof sequencies in time
%     - pss_step:  distance between PSS in I/Q samples
%     - search_window: range around peak to make analysis
%     - correlation_factor: value to determine if there is a peak or not
%     - reample_factor: resample factor to apply
%     - flip: (debug parameter) to flip the signal and zadoof sequencies
%
% Output:
%     - PPM: value of drift
%     - PSS: % of PSS detected
%     - data -> (1) Number of the PSS (iteration) (2) -> Cummulative drift
%     - Y: polyval output


% TODO:
%   - Resampling in the training area

function [PPM, PSS, data, Y, p, PPM2, th, Z_sequence, last, p_loc20] = getDrift (capbuf, sampling_rate, Z, Z_t, PREAMBLE, pss_step, search_window, correlation_factor, resample_factor, flip, degree)
    
    N = 128;
    TRAINING_SAMPLES = pss_step*PREAMBLE;
    
    %TODO: Add upsampling in the training area.
    chunk = capbuf(:);
    
    if (size(chunk,1) < TRAINING_SAMPLES)
        fprintf('[ERROR] Not enought IQ samples, skipping chunk\n');
        [PPM, PSS, data, Y, p, PPM2, th, Z_sequence, last, p_loc20] = deal(NaN);
        return;
    end
    chunk_p = chunk(1:TRAINING_SAMPLES);
    
    p_value = NaN;
    p_loc = NaN;
    Z_sequence = 0;
    max_corr = 0;
    th = NaN;
    % Training step 
    % Check which Zadoff sequence gets higher correlation
    
    for i=1:length(Z_t)
        
        % Size of mask
        mshort = Z_t{i};
        power_mshort=sqrt( sum(abs((mshort)).^2) );

        xcor0 = abs(xcorr(chunk_p,mshort));
        xcor0 = xcor0(size(chunk_p):end);
       
        tmpnorm2 = sqrt(movsum((abs(chunk_p)).^2, [1 N]));                  % power of local signal        
        xcor_final = xcor0(1:size(tmpnorm2,1)) ./ (tmpnorm2*power_mshort);  % normalize the cross-corr
                       
        % compute threshold
        tmp_sorted=sort(xcor_final(1:TRAINING_SAMPLES),'descend');
        th_learned=(tmp_sorted(PREAMBLE)+tmp_sorted(PREAMBLE*2+1))/2;
        th_learned=th_learned*0.70;                        
 
        [p_value_aux, p_loc_aux] = findpeaks(xcor_final,'MinPeakDistance',pss_step-search_window,'MinPeakHeight', th_learned);
        
        if (max(p_value_aux)>max_corr)
            Z_sequence = i;
            p_value = p_value_aux;
            p_loc = p_loc_aux;            
            max_corr = max(p_value_aux);
            
            % compute threshold
            th = th_learned;
        
        end
    end
 
    
    if (size(p_loc,1)) == 0
      fprintf('[ERROR] No find peaks in 1st step, skipping chunk\n');
      [PPM, PSS, data, Y, p, PPM2, th, Z_sequence, last, p_loc20] = deal(NaN);
      return;
    end
    diff_peaks = diff(p_loc);
    diff_peaks = abs(diff_peaks - pss_step);

    if (size(diff_peaks,1) <= 2)
      fprintf('[ERROR] No enough peaks, skipping chunk\n');
      [PPM, PSS, data, Y, p, PPM2, th, Z_sequence, last, p_loc20] = deal(NaN);
      return;
    end

    dist_peaks = find(diff_peaks(end-2:end) > 10 );
    if (size(dist_peaks) ~= 0 )
       fprintf('[WARNING] Some PSS detected are further than %d +- 10 I/Q samples\n', pss_step);      
    end

    last_valid = find (diff_peaks < 10);

    if (size(last_valid,1) == 0)
       fprintf('[ERROR] No valid PSS at the begining, skipping file\n');
       [PPM, PSS, data, Y, p, PPM2, th, Z_sequence, last, p_loc20] = deal(NaN);
       return;
    end
    
    

    SAFE_GUARD=0;

    [p_value20, p_loc20] = getPEAKS(chunk(p_loc(last_valid(end))+SAFE_GUARD:end) , ... 
        pss_step, search_window, th, resample_factor, Z{Z_sequence}, Z_t{Z_sequence}, flip, NaN);

    last = p_loc20(end) + p_loc(last_valid(end));

    [data, Y, p] = analyzeDrift(p_loc20, pss_step, degree);       
    PPM=((Y(end)-Y(1))/((data(end,1)-data(1,1))*pss_step))*1e6;
    
    PPM2 = (p(end-1)/pss_step)*1e6; % Another method (testing purposes)
    
    if (abs(PPM) > 100)
       [PPM, PSS, data, Y, p, PPM2] = deal(NaN);
    else
        PSS = (size(data,1)*100)/((size(chunk(p_loc(last_valid(end))+SAFE_GUARD:end),1)/sampling_rate)*(sampling_rate/pss_step));



    end

    
    

    
    