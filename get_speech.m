function signal = get_speech(spect_mag,spect_phase,f,freq_res,window,stride,use_hamm)
    window_len = uint64(f*window);
    overlap_len = uint64(f*(window-stride));
    signal = [];
    spect_mag = [spect_mag; flipud(spect_mag(2:end-1,:))];
    spect_phase = [spect_phase; flipud(-1*spect_phase(2:end-1,:))];
    
    if use_hamm
        inverse_hamm = 1./hamming(window_len);
    end
    
    for i = 1:size(spect_mag,2)
        fft_seg = spect_mag(:,i) .* exp(1j*spect_phase(:,i));
        signal_segment = real(ifft(fft_seg,freq_res));
        
        if use_hamm
            signal_segment = inverse_hamm.*signal_segment(1:window_len);
        else
            signal_segment = signal_segment(1:window_len);
        end
        
        if isempty(signal)
            signal = signal_segment;
        else
%             disp([length(signal((end-overlap_len)+1:end)), length(signal_segment(1:overlap_len))]);
            overlap_part = (signal((end-overlap_len)+1:end) + signal_segment(1:overlap_len)) / 2;
            non_overlap_part = signal_segment(overlap_len+1:end);
            signal(end-overlap_len+1:end) = overlap_part;
            signal = [signal;non_overlap_part];
        end
    end
end