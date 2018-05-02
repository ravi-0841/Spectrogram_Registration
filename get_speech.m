function signal = get_speech(spect_mag,spect_phase,f,freq_res,window,overlap)
    window_len = f*window;
    overlap_len = f*overlap;
    signal = [];
    spect_mag = [spect_mag; flipud(spect_mag(2:end-1,:))];
    spect_phase = [spect_phase; flipud(-1*spect_phase(2:end-1,:))];
    for i = 1:size(spect_mag,2)
        fft_seg = spect_mag(:,i) .* exp(1j*spect_phase(:,i));
        signal_segment = real(ifft(fft_seg,freq_res));
        signal_segment = signal_segment(1:window_len);
        if isempty(signal)
            signal = signal_segment;
        else
            overlap_part = (signal(end-overlap_len+1:end) + signal_segment(1:overlap_len)) / 2;
            non_overlap_part = signal_segment(overlap_len+1:end);
            
            signal(end-overlap_len+1:end) = overlap_part;
            signal = [signal;non_overlap_part];
        end
    end
end