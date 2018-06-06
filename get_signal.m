function signal = get_signal(mstft, W, S, iterations, wshift, look_ahead)

    if nargin<6;    look_ahead = 3;        end
    if nargin<5;    wshift = 128;          end
    if nargin<4;    iterations = 1000;     end
    if nargin<3;    error;                 end
    
    L=5;
    weights=create_weights(W,S,wshift,L);

    %% Create weights for asymmetric windows
    [W_asym_init,W_asym_full] = build_asymmetric_windows(W.*S,wshift);
    weights_asym_init=create_weights(W_asym_init,S,wshift,L);
    weights_asym_full=create_weights(W_asym_full,S,wshift,L);
    
    online_iterations = iterations;
    online_alpha      = 1;
    online_beta       = 0.1;
    online_gamma      = 1;
    online_thresholds = online_alpha*exp(-online_beta*(0:(online_iterations-1)).^online_gamma);

    X1=online_lws(mstft,weights,weights_asym_init,weights_asym_full,online_thresholds,look_ahead);

    %% run batch LWS
    alpha      = 100;
    beta       = 0.1;
    gamma      = 1;
    thresholds = alpha*exp(-beta*(0:(iterations-1)).^gamma);

    Y=batch_lws(X1,weights,thresholds);
    signal=istft(Y,wshift,S);
end