function [spect_mag,spect_phase] = get_spectrogram(signal,fs,freq_res,window,stride,nmz)
    signal = reshape(signal, [length(signal),1]);
    sample_size = fs*window;
    sample_stride = fs*stride;
    spect_mag = [];
    spect_phase = [];
    i = 1;
    for n_frames=1:ceil((length(signal) - sample_size) / sample_stride + 1)
        try 
            signal_segment = signal(i:i+sample_size-1);
        catch IE
            signal_segment = signal(i:end);
            signal_segment = [signal_segment; zeros(sample_size-length(signal_segment),1)];
        end
        segment = signal_segment;
        i = i + sample_stride;
        fft_segment_mag = abs(fft(segment, freq_res));
        fft_segment_phase = angle(fft(segment, freq_res));
        if iscolumn(fft_segment_mag)
            spect_mag = [spect_mag fft_segment_mag];
            spect_phase = [spect_phase fft_segment_phase];
        else
            spect_mag = [spect_mag fft_segment_mag'];
            spect_phase = [spect_phase fft_segment_phase'];
        end
    end
    spect_mag = spect_mag(1:freq_res/2+1, :);
    spect_phase = spect_phase(1:freq_res/2+1, :);
    
    if nmz==1
        spect_mag = normc(spect_mag);
    end
end