function baseShifter = generateBaseShifterSigmoid(dataStructure)
%   baseShifter = generateBaseShifter(fftl,transitionWidth)

%   linear phase shifter for fractional pitch implementation
%   Designed and coded by Hideki Kawahara
%   24/Feb./2012
%   04/Mar./2012

fftl = (size(dataStructure.spectrogramSTRAIGHT,1)-1)*2;
transitionWidth = dataStructure.transitionWidth;
baseShifter = (0:fftl-1)'/fftl*2*pi;
baseShifter(baseShifter>pi) = baseShifter(baseShifter>pi)-2*pi;
shaper = hanning(round(fftl*transitionWidth/2)*2+1);
shaperBase = -round(fftl*transitionWidth/2):round(fftl*transitionWidth/2);
baseShifter(shaperBase+1+fftl/2) = baseShifter(shaperBase+1+fftl/2).*sqrt(1-shaper);
return;