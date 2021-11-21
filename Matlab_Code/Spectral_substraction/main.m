close all;
clear all;
clc;
pkg load signal;

R=[];
Paco=[];
Y=[];
[X,fs]  =  audioread("fem02SW_c31_bgnpubtk1_Mic07ContactMic.wav");
##[Y,R,Paco] = fw_hpf_vad_aco(X,fs);


[ss,po]=specsubm_mod(X,fs);
audiowrite("out.wav",real(ss),fs);
audiowrite("in.wav",X,fs);

plot(X,'r'); hold on;plot(real(ss));