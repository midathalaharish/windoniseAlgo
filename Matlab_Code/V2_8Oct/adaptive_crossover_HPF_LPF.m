%adaptive_crossover_HPF_LPF


clc; close all; clear;

[x_in, Fs ] = audioread(input);

%% pass input through HPF of fc 80

[b, a]  =  butter(6, 80/Fs/2, 'high');
x_right_acoustic = filter(b, a, x_in(:,1));
x_left_vpu = filter(b, a, x_in(:,2));

%% pass input through BPF for vpu only
band_freqs = [100, 1500];
[b_high, a_high]  =  butter(2, band_freqs(1)/Fs, 'high');
[b_low, a_low]  =  butter(2, band_freqs(2)/Fs, 'low');
b2 = conv(b_high, b_low);
a2 = conv(a_high, a_low);
x_left_vpu_filteed = filter(b2, a2, x_in(:,2));

%% start algorithm

%% declare variables

noverlap  =   256;
window    =   hann(512, 'periodic');
M         =   length(window);
w         =   window(:);
nhop      =   M-noverlap;
nx        =   length(x_left_vpu);
nframes   =   1+fix((nx-M)/nhop);
xframe_l  =   zeros(M,1);
xframe_r  =   zeros(M,1);
exp_avg_r =   0;
NE        =   zeros(M,1);
exp_avg_r_out = zeros(nframes,1);
times     =   (0:length(nframes)-1)*1/Fs;
yout      =   zeros(nframes*nhop+M,2);
times_out =   (0:length(yout)-1)*1/Fs;
filtered_output = zeros(nframes*nhop+M,2);
vad_signal = zeros(M,1);
VAD_signal1 = zeros(nframes*nhop+M,1);
E_out     =   zeros(M,1);
vad_threshold = 2e-5;
E_full    =   zeros(nframes,1);
exp_avg_r_out_dB = zeros(nframes*nhop+M,1);
sum_ne    = zeros(M,1);
Wn        = ones(nframes,1);