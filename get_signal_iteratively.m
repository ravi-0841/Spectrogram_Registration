function est_signal = get_signal_iteratively(mag, N, wshift, max_iter)
    W=sqrt((0.5-0.5*cos(2*pi*(0:(N-1))'/(N)))/(N/wshift)*2);
    est_signal = istft(mag,wshift,W);
    X = stft(est_signal,N,wshift,W);
    for i = 1:max_iter
        est_signal = istft(mag.*exp(1j*angle(X)), wshift, W);
        X = stft(est_signal,N,wshift,W);
    end
end