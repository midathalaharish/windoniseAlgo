function [x,fs] = se_input(wav_vpu,wav_mic)
% [x,fs] = se_input(wav_vpu,wav_mic)
% [x,fs] = se_input

if nargin==0
    wav_vpu = uigetfullfile('*.wav','Select Contact Mic Wav-File to Open');    
    wav_mic = uigetfullfile('*.wav','Select Acoustic Mic Wav-File to Open',fileparts(wav_vpu));    
end

info_vpu = audioinfo(wav_vpu);
info_mic = audioinfo(wav_mic);

if info_vpu.SampleRate ~= info_vpu.SampleRate
    error('The wav-files must have the same sample rate.')
end

if (info_vpu.NumChannels>1) || (info_mic.NumChannels>1)
    error('The wav-files must be mono (1 channel).')
end

[Y_vpu,fs] = audioread(wav_vpu);
Y_mic = audioread(wav_mic);

n = min(length(Y_vpu),length(Y_mic));
x = [Y_vpu(1:n) Y_mic(1:n)];

end