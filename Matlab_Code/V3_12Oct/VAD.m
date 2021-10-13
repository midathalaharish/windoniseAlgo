%%% WOLA VAD

function [VADOut, E_out] = VAD(xsample, E, tA, M, nhop, vad_threshold, Fs)
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
alpha = 1/(Fs*tA);
const = 0.6;
hang_time = 5e-2;
hang_samples  = round(hang_time * Fs);
count =0;
Start_skip_samples =0;
%release = 0.7;
%attack = 0.001;
%    NoiseMargin=3;
%    Hangover=8;
%    NoiseCounter=0;
%    SpeechCounter=0;
p = E;                  % start with nhop(256) and read the last updated value
for k = 1:M
% Compute sample energy
if(Start_skip_samples == 1)
   VADOut(k) = 1;
   count =count +1;
##    p = (1-alpha)*p + alpha*xsample(k).^2; ##
     E_out(k) = p;
   if(count == hang_samples) || (k == M)
      Start_skip_samples=0;
   endif
   continue
endif
  p = (1-alpha)*p + alpha*xsample(k).^2;
  E_out(k) = p;
  
  if (E >= vad_threshold)
    VADOut(k) = 1;
    Start_skip_samples=1;
    
  else
    VADOut(k) = 0;
    Start_skip_samples=0;
  endif 


if(0)
for k = 1:M
% Compute sample energy
  p = (1-alpha)*p + alpha*xsample(k).^2;
  E_out(k) = p;
  if (E >= vad_threshold || ((VAD_delay > 0.7))) 
    VADOut(k) = 1;
    VAD_delay1 = 1;
  else
    VADOut(k) = 0;
    VAD_delay1 = VAD_delay*const;
    if(VAD_delay1 < 0)
      VAD_delay1 = 0;
     endif 
  end
end
endif
end