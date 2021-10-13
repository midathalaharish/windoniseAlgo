function [VADOut, E_out] = VAD(xsample, E, tA, M, nhop, vad_threshold)
%**************************************************************************
% Voice activity detector (VAD)
%
% Implements a VAD that is based on signal to noise ratio (SNR) estimates.
% Subband SNR estimates are thresholded towards a preset threshold and a
% majority decision is taken to form the detection. The VAD function should
% be called once for every new input sample.
%
% Input parameters:
%  - xsample (scalar) the input signal sample value
%  - E (Kx1 vector) energy estimate
%  - Eflr (Kx1 vector) energy estimate
%  - tA (scalar) time constant A
%  - tB (scalar) time constant B
%  - tBF (scalar) time constant BF
%  - SNR_THS (scalar) signal to noise ratio threshold
%  - B (Kx2 matrix) VAD analysis filter bank FIR coefficients
%  - A (Kx3 matrix) VAD analysis filter bank IIR coefficients
%  - n (scalar) sample index

%
% Output parameters:
%  - VADOut (scalar) VAD detector output (1=speech, 0=noise)
%  - E (Kx1 vectir) updated energy estimate
%  - Eflr (Kx1 vectir) updated energy estimate

%**************************************************************************

p = E;                        
for k = 1:M

% Compute sample energy
  p = (1-tA)*p + tA*xsample(k).^2;
  E_out(k) = p;
  if (E >= vad_threshold)
    VADOut(k) = 1;
  else
  VADOut(k) = 0;
  end
end