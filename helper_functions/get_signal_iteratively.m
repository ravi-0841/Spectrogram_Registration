function est_signal = get_signal_iteratively(mag, N, wshift, W, max_iter)
    
    if nargin<5;    max_iter = 1000;    end
    
    est_signal = istft(mag,wshift,W);
    X = stft(abs(est_signal),N,wshift,W);
    for i = 1:max_iter
        est_signal = istft(mag.*exp(1j*angle(X)), wshift, W);
        X = stft(est_signal,N,wshift,W);
    end
end