function responseInFreqDomain = minimumPhaseResponse(specSlice)
%   responseInFreqDomain = minimumPhaseResponse(specSlice)

%   minimum phase resoponse through complex cepstrum
%   codec by Hideki Kawahara
%   24/Feb./2012
%   04/Mar./2012 revised to generalize
%   01/Aug./2017 bug fix (safeguard)

%responseInFreqDomain = [];
fftl = 2*(size(specSlice,1)-1);
specSlice(specSlice == 0) = min(specSlice(specSlice > 0)); % safe guard
doubleSpectrum = [specSlice;specSlice(end-1:-1:2)];
complexCepstrum = ifft(log(doubleSpectrum)/2);
complexCepstrum(fftl/2+1:end) = 0;
complexCepstrum(2:fftl/2) = complexCepstrum(2:fftl/2)*2;
responseInFreqDomain = exp(fft(complexCepstrum));
return;