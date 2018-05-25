function [i_mag, i_phase] = interpolate_stft_from_dtw(mag,phase,p)
    i_mag = zeros(size(mag,1),length(p));
    i_phase = zeros(size(phase,1),length(p));
    
    i_mag(:,1) = mag(:,1);
    i_phase(:,1) = phase(:,1);
    
    counter = 2;
    
    uniq_vals = unique(p);
    for value = 1:length(uniq_vals)
        repeats = find(p==uniq_vals(value));
        if length(repeats)>1
            start_val = uniq_vals(value);
            end_val = uniq_vals(min([value+1, length(uniq_vals)]));
            
            i_mag(:,counter) = mag(:,start_val);
            i_phase(:,counter) = mag(:,start_val);
            counter = counter  + 1;
            
            delta_mag = (mag(:,end_val) - mag(:,start_val)) / length(repeats);
            delta_phase = (phase(:,end_val) - phase(:,start_val)) / length(repeats);
            
            for r = 1:length(repeats)-1
                disp([start_val, end_val]);
                i_mag(:,counter) = mag(:,start_val) + r*delta_mag;
                i_phase(:,counter) = phase(:,start_val) + r*delta_phase;
%                 z = mvnrnd(phase(:,end_val), eye(257,257));
%                 i_phase(:,counter) = z';
                counter = counter + 1;
            end
        else
            i_mag(:,counter) = mag(:,uniq_vals(value));
            i_phase(:,counter) = phase(:,uniq_vals(value));
            counter = counter + 1;
        end
    end
%     figure()
%     subplot(211), plot(i_phase(10,:)), title('Phase')
%     subplot(212), plot(i_mag(10,:)), title('Magnitude')
end