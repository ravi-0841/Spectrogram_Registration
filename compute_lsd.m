function Dls = compute_lsd(source,target)
    source_fft = abs(fft(source, 512)).^2;
    target_fft = abs(fft(target, 512)).^2;
    Dls = sqrt(1/length(source_fft) * sum(10*log10(source_fft./target_fft)));
end