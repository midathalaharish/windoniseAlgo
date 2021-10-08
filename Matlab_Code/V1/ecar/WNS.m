%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:  WNS_ALgorithm.m
%
% Description: Suppres wind noise by dual mic RMS_Diff WND method. 
%              Once detected wind noise, will decrease the gain in that frame. 
%              Haven't added the hang over and smooth parts!!!
%
% Based on:    WND at same folder.
%
% Input file:  mic6_mic7_006A_fan3ms_speech.wav a 6 minute long file with
%              two channel data, ch1 is acoustic and ch2 is contact mic 
%
% Date:        8/18/2021, initial
% Modified:    8/28/2021, suppress with new OTOG2 data from Gongqiang
% Modified:    9/03/2021, coninue with OTOG2 data
% Note, modify input, output, MSC and RMSD threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean up
clc;
clear; 
close all;

%-----------------------------------------------------------------------
% Apply WND to real contact mic data  
%
% From Lalin 6 minutes long data speech + wind 3m/s. 
% Ch1 is from acoustic mic, Ch2 is from contact mic
%-----------------------------------------------------------------------

%-------------------------------------------------
% Lalin's long file
%-------------------------------------------------
% set common folder for wav files
%folderName = 'C:\Projects\AfterAug2\DataFromLalin\';
% set file names
%fileName = strcat(folderName, 'mic6_mic7_006A_fan3ms_speech.wav'); % Acoustic + VPU
%[x_in, Fs] = audioread(fileName);
%inputInfo  = audioinfo(fileName)         % two channels, 16000s/s, 16bit/s
% get contact mic VPU signal and acoustic mic signal, add some gain!!!
%x_left_vpu = x_in(:,2)*20;                %10
%x_right_acoustic = x_in(:,1)*20;          %10


%------------------------------------------------------
% 9/3/2021, continue OTOG2 data
%------------------------------------------------------
% set common folder for wav files
folderName = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\';

%-----------------------------------------------------------------------
% OTOG2_case01: speech_wind_3mps_FanB
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');

%-----------------------------------------------------------------------
% OTOG2_case02: speech_wind_3mps_TunnelA
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');

%-----------------------------------------------------------------------
% OTOG2_case03: speech_wind_5mps_FanA
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');

%-----------------------------------------------------------------------
% OTOG2_case04: speech_wind_5mps_TunnelB
%-----------------------------------------------------------------------
% set file name
fileName = strcat(folderName, 'OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');


% get data
[x_in, Fs] = audioread(fileName);
audioinfo(fileName)                    % two channels, 48000s/s, 16bit/s
% get contact mic VPU signal and acoustic mic signal, add some gain!!!
x_left_vpu = x_in(:,2);
x_right_acoustic = x_in(:,1);


%-------------------------------------------------------------------------
% start WNS
%-------------------------------------------------------------------------

% get time range values in seconds
times = (0:length(x_left_vpu)-1)*1/Fs;          % just for time domain plot

% set up signal analysis parameters
noverlap=256;  % number of samples the sections of A overlap. If negative, 
               % -NOVERLAP is the "hop size". (The overlap is the window 
               % length minus the hop size.) NOVERLAP < M
%window=hamming(512);  % length of M window function, in ZERO-PHASE form,
                      % to each frame of A. M < NFFT, NFFT-M zero pad
window = hann(512, 'periodic');                      
nfft=2048;     % number point of FFT, NFFT. 1024 or 2048? will test

alpha = 0.91;  % coherence smooth factor 0.8

% This script splits the input signal into overlapping segments of
% length M, windows each segment with the length M WINDOW vector.

% check window
M = length(window);
Modd = mod(M,2);   % 0 if M even, 1 if odd
Mo2 = (M-Modd)/2;
w = window(:);     % Make sure it's a column

% set hop number, assume 0<noverlap<M
nhop = M-noverlap;

% set input data
nx = length(x_left_vpu);            % assume length(x) > M
% another way to calculate frame numbers
nframes = 1+fix((nx-M)/nhop);

X_left  = zeros(nfft,1);             % allocate output spectrum in advance
X_right = zeros(nfft,1);

SX_left  = zeros(nfft,nframes);      % allocate output spectrogram in advance
SX_right = zeros(nfft,nframes);

MSC_SX  = zeros(nfft,nframes);
SSC_SX  = zeros(nfft,nframes);       % allocate VPU SSC spectrogram in advance
RMSD_SX = zeros(nfft,nframes);       % allocate RMSD spectrogram in advance

zp = zeros(nfft-M,1);         % zero-padding for each FFT
xframe_l = zeros(M,1);        % each frame samples
xframe_r = zeros(M,1);

% smoothed auto, cross power spectrum and MSC
XC_11 = zeros(nfft,1);
XC_22 = zeros(nfft,1);
XC_12 = zeros(nfft,1);
MSC = zeros(nfft,1);

% WND, coherence method
wind_decision = zeros(nx,1); 
detect_result = zeros(nx,1);
th_msc = 0.79;
MSC_rate = zeros(nframes,1);

% WND, SSC method
SSC_decision_1 = zeros(nx,1);
SSC_decision_2 = zeros(nx,1);
SSC_rate = zeros(nframes,1);

% RMS diff
rms_diff = zeros(M,1);
%%%rms_diff = 0;
RMSD_decision = zeros(nx,1);
RMSD_rate = zeros(nframes,1);

AllT_rate = zeros(nframes,1);

% processed/enhanced audio data
yout = zeros((nframes-1)*nhop+M, 2);

times_out = (0:length(yout)-1)*1/Fs;

counter = 0;

% frame processing
for m=0:nframes-1
  xframe_l = x_left_vpu(1+m*nhop : M+m*nhop);
  xframe_r = x_right_acoustic(1+m*nhop : M+m*nhop);
  
  % windowing and fft
  xw_l = w .* xframe_l;                         % Apply window
  xw_r = w .* xframe_r;
  xwzp_l = [xw_l(Mo2+1:M);zp;xw_l(1:Mo2)];      % in zero-phase form!!
  xwzp_r = [xw_r(Mo2+1:M);zp;xw_r(1:Mo2)];
  %xwzp_l = [xw_l;zp];        % be careful this for frequency domain processing
  %xwzp_r = [xw_r;zp];
  X_left  = fft(xwzp_l);
  X_right = fft(xwzp_r);
  
  % spectrogram of input signal
  SX_left(:,m+1) = X_left;
  SX_right(:,m+1) = X_right;
  
  % coherence  
  XC_11 = alpha*XC_11 + (1-alpha)*(X_left.*conj(X_left));
  XC_22 = alpha*XC_22 + (1-alpha)*(X_right.*conj(X_right));
  XC_12 = alpha*XC_12 + (1-alpha)*(X_left.*conj(X_right));
  % MSC
  MSC = (XC_12.*conj(XC_12)) ./ (XC_11.*XC_22);
  
  MSC_SX(:,m+1) = MSC;
  
  % average of lower frequency 0 - 500Hz, 500/(Fs/nfft) -> 21.xx -> 22
  MSC_ave = sum(MSC(1:44))/44;
  wind_decision(m*nhop+(1:nhop),1) = 1 - MSC_ave;
  % set flags
  if ((1 - MSC_ave) > th_msc)
      MSC_rate(m+1,1) = 1; 
%%%      detect_result((m)*nhop+(1:nhop),1) = 1;
  else
%%%      detect_result((m)*nhop+(1:nhop),1) = 0.5;
      MSC_rate(m+1,1) = 0;
  end 
  
  % SSC
  SSC_n_1 = 0;
  SSC_d_1 = 0;
  SSC_n_2 = 0;
  SSC_d_2 = 0;
  for mm = 1:20
      SSC_n_1 = SSC_n_1 + XC_11(mm)*mm*Fs/nfft;
      SSC_d_1 = SSC_d_1 + XC_11(mm);
      SSC_n_2 = SSC_n_2 + XC_22(mm)*mm*Fs/nfft;
      SSC_d_2 = SSC_d_2 + XC_22(mm);
  end
  SSC_1      = SSC_n_1 / SSC_d_1;
  SSC_1_flag = ((20*Fs/nfft) - SSC_1) / (20*Fs/nfft);
  SSC_decision_1((m)*nhop+(1:nhop),1) = SSC_1_flag;
  SSC_2      = SSC_n_2 / SSC_d_2;
  SSC_2_flag = ((20*Fs/nfft) - SSC_2) / (20*Fs/nfft);
  SSC_decision_2((m)*nhop+(1:nhop),1) = SSC_2_flag;
  
  SSC_SX(:,m+1) = ones(nfft,1)*SSC_1_flag;

  % use VPU only
  if (SSC_1_flag > 0.55)     % 0.7, 0.4
      SSC_rate(m+1,1) = 1;
%%%      detect_result((m)*nhop+(1:nhop),1) = 1;    
  else
      SSC_rate(m+1,1) = 0;
%%%      detect_result((m)*nhop+(1:nhop),1) = 0.5;
  end

  % RMS diff
  rms_left    = rms(xframe_l);
  rms_dB_left = 20*log10(rms_left+1.0e-12);     % +1.0e-12 to avoid log(0)!!!
  rms_right    = rms(xframe_r);
  rms_dB_right = 20*log10(rms_right+1.0e-12);

  rms_diff = 0.95*rms_diff + (1-0.95)*(rms_dB_right - rms_dB_left);    %0.98, 0.85 for Long file
  RMSD_decision((m)*nhop+(1:nhop),1) = mean(rms_diff);
  
  RMSD_SX(:,m+1) = ones(nfft,1)*mean(rms_diff);

  % use RMS diff only
  if (mean(rms_diff) > 14.2)      %18.5 for Sigma data, 7.2 for Sage data, 14.5 for Long file
      RMSD_rate(m+1,1) = 1;
%%%      detect_result((m)*nhop+(1:nhop),1) = 1;
  else
      RMSD_rate(m+1,1) = 0;
%%%      detect_result((m)*nhop+(1:nhop),1) = 0.5;
  end
  
  % MSC, SSC_VPU and RMS Diff
  
  %%%AA = ((1 - MSC_ave) > th_msc);    % for Long file
  AA = ((1 - MSC_ave) > 0.88);      % 9/3/2021, 0.92 for OTOG2 case01, 0.90 for case2
  BB = (SSC_1_flag > 0.4);
  %%%CC = (mean(rms_diff) > 14.5);     % 18.5 for Sigma data, 7.2 for Sage data, 14.5 for Long file
  
  CC = (mean(rms_diff) > 16.5);
  %if ~(CC)     % for Long file
  %if ~(AA)     % for first 3 cases
  if ~(CC)
      detect_result((m)*nhop+(1:nhop),1) = 0.5;
      AllT_rate(m+1,1) = 0;
      
      frameGain = 1.0;
      counter = 0;
  %elseif AA || BB || (CC)   % for Long file
  else   
      detect_result((m)*nhop+(1:nhop),1) = 1.0;
      AllT_rate(m+1,1) = 1;
      
      counter = counter + 1;
      if (counter >= 10)
         frameGain = 0.15;    %0.1
      end
  end

  % processed output
  %y_left = (xframe_l / 20) * frameGain;            % for Long file /20
  y_left = (xframe_l / 2.0) * frameGain;
  yout(m*nhop+(1:M), 2) = yout(m*nhop+(1:M), 2) + y_left;   % Ch2 to match with input
  %y_right = (xframe_r / 20) * frameGain;           % for Long file /20
  y_right = (xframe_r / 2.0) * frameGain;
  yout(m*nhop+(1:M), 1) = yout(m*nhop+(1:M), 1) + y_right;
   
end

%------------------------------------------------------
% write output data to a file
%------------------------------------------------------
% Lalin's data
%outputFileName = strcat(folderName, 'Output\wns_mic6_mic7_006A_fan3ms_speech.wav');
%audiowrite(outputFileName, yout, Fs);
%audioinfo(outputFileName)

% OTOG2_case01:
%outputFileName = strcat(folderName, 'WindSuppressed\wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
% OTOG2_case02:
%outputFileName = strcat(folderName, 'WindSuppressed\wns_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
% OTOG2_case03:
%outputFileName = strcat(folderName, 'WindSuppressed\wns_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
% OTOG2_case04:
outputFileName = strcat(folderName, 'WindSuppressed\wns_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');


% save data
audiowrite(outputFileName, yout, Fs);
audioinfo(outputFileName)



  % spectrograms
  t = (M/2:nhop:M/2+(nframes-1)*nhop)/Fs;    % time in seconds
  f = 0.001*(0:nfft/2)*Fs/nfft;              % frequency
  Xdb_l = 20*log10(abs(SX_left));
  Xdb_r = 20*log10(abs(SX_right));
  Xdb_left = Xdb_l(1:nfft/2+1, :);
  Xdb_right = Xdb_r(1:nfft/2+1, :);
  
  MSC_SX_half = MSC_SX(1:nfft/2+1, :);
  SSC_SX_half = SSC_SX(1:nfft/2+1, :);
  RMSD_SX_half = RMSD_SX(1:nfft/2+1, :);
  
  % set up figure
  position = [50, 50, 900, 600];
  set(gcf, 'position', position)

  figure (1);
  subplot(2,1,1);
  plot(times, x_right_acoustic', 'b');
  set(gca, 'FontName', 'Arial', 'FontSize', 12);   %'linewidth', 2
  hold on;
  grid on;
  plot(times, x_left_vpu', 'y');
  plot(times, detect_result, 'g--', 'linewidth', 1.0);
  xlabel('Time (Second)');
  ylabel('Amplitude');
  title('Time domain voice and wind signals');
  %%%axis([0 36 -1.1 1.1]);
  axis([0 20 -1.1 1.1]);
  %%%axis([0 6 -1.1 1.1]);
  legend('Acoustic microphone', 'Contact microphone');

  subplot(2,1,2);
  % get rid of first 380 samples
  plot(times, wind_decision*20, 'r', 'linewidth', 0.5);
%  plot(times, SSC_decision_1, 'm', 'linewidth', 0.5);
%  plot(times, RMSD_decision, 'm', 'linewidth', 0.5);

  set(gca, 'FontName', 'Arial', 'FontSize', 12);
  grid on;
  hold on;
  
  plot(times, SSC_decision_1*20, 'c', 'linewidth', 0.5);
  %plot(times, SSC_decision_2, 'c', 'linewidth', 0.5);
  
  plot(times, RMSD_decision, 'g', 'linewidth', 0.5);
  
    
  xlabel('Time (Second)');
  ylabel('Amplitude');
  title('Wind detection indicator');
  %axis([0 36 -0.1 1.1]);
  %%%xlim([0 36]);
  xlim([0 20]);
  %%%xlim([0 6]);
  legend('MSC', 'SSC VPU', 'RMS Diff');
  %legend('RMS diff');
  %legend('MSC');
  %legend('SSC VPU');
 
  %figure (2);
  %plot(times, x_right_acoustic', 'b');
  %set(gca, 'FontName', 'Arial', 'FontSize', 12);   %'linewidth', 2
  %hold on;
  %grid on;
  %plot(times, x_left_vpu', 'g');
  %xlabel('Time (Second)');
  %ylabel('Amplitude');
  %title('Time domain voice and wind signals');
  %%%axis([0 36 -1.1 1.1]);
  %axis([0 6 -1.1 1.1]);
  %legend('Acoustic microphone', 'Contact microphone');
  
  figure(3);
  subplot(5,1,1);
  imagesc(t,f,Xdb_left);
  % grid;
  axis('xy');
  colormap(jet);      % low frequency content of the first portion of signal 
                      % is displayed in the lower left corner of the axes.
  colorbar;
  xlabel('Time (Second)');
  ylabel('Frequency (KHz)');
  title('Spectrogram of contact mic signal');
  
  subplot(5,1,2);
  imagesc(t,f,Xdb_right);
  % grid;
  axis('xy');
  colormap(jet);      % low frequency content of the first portion of signal 
                      % is displayed in the lower left corner of the axes.
  colorbar;
  xlabel('Time (Second)');
  ylabel('Frequency (KHz)');
  title('Spectrogram of acoustic mic signal');
  
  subplot(5,1,3);
  imagesc(t,f,MSC_SX_half);
  % grid;
  axis('xy');
  colormap(jet);      % low frequency content of the first portion of signal 
                      % is displayed in the lower left corner of the axes.
  colorbar;
  xlabel('Time (Second)');
  ylabel('Frequency (KHz)');
  title('Spectrogram of MSC');
  
  subplot(5,1,4);
  imagesc(t,f,SSC_SX_half);
  % grid;
  axis('xy');
  colormap(jet);      % low frequency content of the first portion of signal 
                      % is displayed in the lower left corner of the axes.
  colorbar;
  xlabel('Time (Second)');
  ylabel('Frequency (KHz)');
  title('Spectrogram of VPU SSC');
  
  subplot(5,1,5);
  imagesc(t,f,RMSD_SX_half);
  % grid;
  axis('xy');
  colormap(jet);      % low frequency content of the first portion of signal 
                      % is displayed in the lower left corner of the axes.
  colorbar;
  xlabel('Time (Second)');
  ylabel('Frequency (KHz)');
  title('Spectrogram of RMS Diff');
