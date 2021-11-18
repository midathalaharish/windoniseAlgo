%%%%% ADAPTIVE CROSSOVER FREQUENCIES OF HPF-LPF
clc;
close all;
clear all;
pkg load signal
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Read input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xin, Fs] = audioread ("OTOG2_Speech_silent_wind3_wind5_FanB_mic7.wav");


[b,a]            = butter(2, 80/(Fs/2), 'high');
x_right_acoustic = filter(b,a, xin(:,1));
x_left_vpu       = filter(b,a, xin(:,1)); ##change to 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    pre treat vpu signal with BPF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band_freqs = [80,1500];
[b2, a2] = butter(1, [band_freqs(1)/(Fs/2),band_freqs(2)/(Fs/2)], 'bandpass');
x_vpu_filtered = filter(b2, a2, xin(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    variables declaration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nooverlap = 256;
Win = hann(512, 'periodic');
Win_len = length(Win);

nhop  = Win_len;
nx = length(x_left_vpu);
nframes = 1+fix((nx-Win_len)/nhop);
xframe_l = zeros(Win_len,1);
xframe_r = zeros(Win_len,1);

VAD_signal = zeros(Win_len,1);
E_out = zeros(Win_len,1);
vad_threshold = 2e-5;

yout = zeros(nframes*nhop+Win_len,1);

Xc = [80,315,630,1250,2500,5000];

tA = 1/Fs*0.01;
E =0;
exp_avg_r =0;

for m = 0:nframes - 1
     xframe_l           = x_left_vpu(m*nhop+(1:Win_len));
     CCC = lpc(xframe_l,length(xframe_l)-1);
     yout(m*nhop+(1:Win_len), 1)  = yout(m*nhop+(1:Win_len), 1) + CCC';

endfor
     out1 = yout./max(yout);
     plot(x_left_vpu,'b')
     hold on;
     plot(out1,'r');
     grid on;

