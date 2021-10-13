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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    pre treat input with HPF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

nhop  = Win_len-nooverlap;
nx = length(x_left_vpu);
nframes = 1+fix((nx-Win_len)/nhop);
xframe_l = zeros(Win_len,1);
xframe_r = zeros(Win_len,1);

VAD_signal = zeros(Win_len,1);
E_out = zeros(Win_len,1);
vad_threshold = 2e-5;

yout = zeros(nframes*nhop+Win_len,2);

Xc = [80,315,630,1250,2500,5000];

tA = 1/Fs*0.01;
E =0;
exp_avg_r =0;

for m = 0:nframes - 1
 
   xframe_l           = x_left_vpu(m*nhop+(1:Win_len));
   xframe_l_filtered  = x_vpu_filtered(m*nhop+(1:Win_len));
   xframe_r           = x_right_acoustic(m*nhop+(1:Win_len));
    startframeAt = m*nhop

  %--------------------------------------------------------------------------
  %                     VAD
  %--------------------------------------------------------------------------

    [VAD_signal, E_out] = VAD(xframe_l_filtered, E, tA, Win_len, nhop, vad_threshold);
    E = E_out(nhop);      % update E with previous average out value
      
  %--------------------------------------------------------------------------
  %           Adaptive crossover frequencies
  %--------------------------------------------------------------------------
  
   [x_com, xw_l, xw_r, NE,Wn]  =  adaptive_crossover(xframe_l, xframe_r, Xc, Fs, exp_avg_r, VAD_signal, Win_len, Win);
    exp_avg_r = NE(nhop);
##    exp_avg_r_out(m*nhop+(1:nhop), 1)     = NE(1:nhop);
##    exp_avg_r_out_db = 10*log10(exp_avg_r_out+1.0e-12);
    yout(m*nhop+(1:Win_len), 1)   = yout(m*nhop+(1:Win_len), 1) + x_com;
    yout(m*nhop+(1:Win_len), 2)   = yout(m*nhop+(1:Win_len), 2) + xw_r;
    Wn(m+1,1) = Wn;

end
toc