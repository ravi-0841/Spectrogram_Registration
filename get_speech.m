function signal = get_speech(spect_mag,spect_phase,f,window,overlap)
    window_len = f*window;
    overlap_len = f*overlap;
    signal = [];
    for i = 1:size(spect_mag,2)
        fft_seg = spect_mag(:,i) .* exp(1j*spect_phase(:,i));
        signal_segment = ifft(fft_seg,length(fft_seg));
        if isempty(signal)
            signal = signal_segment;
        else
            length((window_len-overlap_len)*(i-1)+1:(i-1)*window_len)
            signal((window_len-overlap_len)*(i-1)+1:(i-1)*window_len) = ...
                (signal((window_len-overlap_len)*(i-1)+1:(i-1)*window_len) + signal(1:overlap_len)) / 2;
            signal((i-1)*window_len+1:i*window_len-overlap_len) = signal(overlap_len+1:end);
        end
    end
end