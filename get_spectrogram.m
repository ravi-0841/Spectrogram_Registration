function [spect_mag,spect_phase] = get_spectrogram(signal,fs,window_size,overlap)
    sample_size = fs*window_size;
    sample_stride = fs*overlap;
    spect_mag = [];
    spect_phase = [];
    for i=1:sample_size-sample_stride:length(signal)-sample_size
        segment = signal(i:i+sample_size);
        fft_segment_mag = abs(fft(segment, sample_size));
        fft_segment_phase = angle(fft(segment, sample_size));
        if iscolumn(fft_segment_mag)
            spect_mag = [spect_mag fft_segment_mag];
            spect_phase = [spect_phase fft_segment_phase];
        else
            spect_mag = [spect_mag fft_segment_mag'];
            spect_phase = [spect_phase fft_segment_phase'];
        end
    end
end