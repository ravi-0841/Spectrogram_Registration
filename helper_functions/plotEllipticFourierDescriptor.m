function plotEllipticFourierDescriptor(a,b,c,d,T,s)
% Plots the elliptic fourier descriptor
% Bibliography:
% F.P. Kuhl and C.R. Giadina, Elliptic Fourier features of a closed
% contour, Comput. Vis. Graphics Image Process. 18, 236-258 (1982)
% W. Dahmen and A. Reusken, Numerik fuer Ingenieure und 
% Naturwissenschaftler, 2. korrigerte Auflage, Springer Verlag, Aachen 2007
% A. Angermann, M. Beuschel, M. Rau and U. Wohlfarth, MATLAB - Simulink - 
% Stateflow, 7. Auflage, Oldenbourg Verlag, M?nchen 2011
% Remark: 
% - We don't use a fast fourier transform. Hence, the calculation will take
%   O(N)
% - There are no parameter checks
%
% Copyright (c) 2013 by Tobias Pohlen (tobias.pohlen@gmail.com)
% www.geekstack.net
    
    % Receive the fourier functions
    xt = ellipticFourierFunction(a,b,T);
    yt = ellipticFourierFunction(c,d,T);
    
    k = 1;
    for it = 0:s:T
        pl(k) = xt(it) + i*yt(it);
        k = k + 1;
    end
    
    plot(pl);
    axis equal;
end