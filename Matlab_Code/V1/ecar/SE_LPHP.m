%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Speech enhancement by HP and LP filter pair.
% Reference paper: Combination of Bone-Conducted Speech with Air-Conducted
%                  Speech Changing Cut-Off Frequency (downloaded)
%                  Distortion Measures for Speech Processing
%
% Test files:  iTP02 - rec_alexa_pizza.wav and 
%              iTP02_Music_rec_alexa_pizza.wav are from Sigma
%
% Date:        4/19/2021, initial, based on VAD_Algorithm_final.m
% Modified:    4/21/2021, add smooth filter for spectrum and distortion 
%                         measure: Log Spectral Deviation.
%              5/4/2021, implement LP model based SE method
%
% Notes:       5/5/2021, need more fine tune and improvement. Check with
%                        the paper and other description. More test and
%                        comparisons
% Modified:    9/4/2021, test with OTOG2 data
% Note:                  cutoff frequency is a key tuning parameter. Lower
%                        it from 8k to 5k to get more signal from acoustic
%                        mic. No need to apply HPF for acoustic mic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean up
clc;
clear; 
pkg load signal
%close all;

%------------------------------------------------------
% 9/4/2021
%------------------------------------------------------
%-------------------------------------------------------------------------
% load input data
%-------------------------------------------------------------------------
%folderName = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\';
%outputFolderName = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\EnhancedData\';

%folderName_wns = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\WindSuppressed\';
% PI
folderName_wns = 'C:\Users\00317859\Desktop\OTO\copy of 2021-09-15 From Lin - sandbox\inputFiles\';
%outputFolderName_ens = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\WindSuppressed\EnhancedData\';
% PI
outputFolderName_ens = 'C:\Users\00317859\Desktop\OTO\copy of 2021-09-15 From Lin - sandbox\outputFiles\';


% load recorded data, a two-channel wave file sampling @ 48 KHz, first, left
% channel is acousitc mic signal, second, right channel is contact mic
% get contact mic VPU signal and acoustic mic signal
% data from Sigma, Note, the order is different
%[x_in, Fs] = audioread('iTP02_Music_rec_alexa_pizza.wav');
%[x_in, Fs] = audioread('iTP02_rec_alexa_pizza.wav');
%x_left_vpu = x_in(:,1);
%x_right_acoustic = x_in(:,2);


%--------------------------------------------------
% 9/4/2021 for new OTOG2 data
%--------------------------------------------------
%-----------------------------------------------------------------------
% OTOG2_case01: speech_wind_3mps_FanB, the input file is created by
%               CreateTwoChannelData.m
%               ch-1: acoustic mic; ch-2: contact mic
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
%fileName_wns = strcat(folderName_wns, 'wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
% PI
%inputFile = 'Stereo_Claes _OTO6.0_2m_cables Mic15-Knowles.wav';
inputFile = 'OTOG2_Speech_silent_wind3_wind5_FanB_mic7.wav';
##fileName_wns = strcat(folderName_wns, inputFile);
fileName_wns= inputFile;

%-----------------------------------------------------------------------
% OTOG2_case02: speech_wind_3mps_TunnelA
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
%fileName_wns = strcat(folderName_wns, 'wns_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');

%-----------------------------------------------------------------------
% OTOG2_case03: speech_wind_5mps_FanA
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
%fileName_wns = strcat(folderName_wns, 'wns_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');

%-----------------------------------------------------------------------
% OTOG2_case04: speech_wind_5mps_TunnelB
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
%fileName_wns = strcat(folderName_wns, 'wns_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');

%-----------------------------------------------
% get data from the file without WNS
%-----------------------------------------------
%[x_in, Fs]  = audioread(fileName);
%inputInfo   = audioinfo(fileName)

%-----------------------------------------------
% get data from the file with WNS
%-----------------------------------------------
[x_in, Fs]  = audioread(fileName_wns);
audioinfo(fileName_wns)

%-----------------------------------------------
% add HPF for acoustic mic
%-----------------------------------------------
[b, a] = butter(6, 300/(48e3/2), 'high');     %9/4/2021
x_right_acoustic = filter(b, a, x_in(:,1));   %9/4/2021

%--------------------------------------------------
% get VPU, acoustic mic data
%--------------------------------------------------
%x_right_acoustic = x_in(:,1);
x_left_vpu = x_in(:,1);


%--------------------------------------------------------------------------
% start algorithm
%--------------------------------------------------------------------------
% get time range values in seconds
times = (0:length(x_left_vpu)-1)*1/Fs;          % just for time domain plot

% set up signal analysis parameters
noverlap=256;  % number of samples the sections of A overlap. If negative, 
               % -NOVERLAP is the "hop size". (The overlap is the window 
               % length minus the hop size.) NOVERLAP < M
%window=hamming(512);  % length of M window function, in ZERO-PHASE form,
                      % to each frame of A. M < NFFT, NFFT-M zero pad
window = hann(512, 'periodic');                      
nfft=2048;     % number point of FFT, NFFT.

%alpha = 0.91;  % smooth factor for coherence
alpha = 0.1;  %0.5;

% This script splits the input signal into overlapping segments of
% length M, windows each segment with the length M WINDOW vector, in
% zero-phase form.

% check window
M = length(window);
Modd = mod(M,2);   % 0 if M even, 1 if odd
Mo2 = (M-Modd)/2;
w = window(:);     % Make sure it's a column

% set hop number, assume 0<noverlap<M
nhop = M-noverlap;

% set input data
nx = length(x_left_vpu);      % assume length(x) > M
% another way to calculate frame numbers
nframes = 1+fix((nx-M)/nhop);

zp = zeros(nfft-M,1);         % zero-padding for each FFT
xframe_l = zeros(M,1);        % framed input samples
xframe_r = zeros(M,1);

X_left  = zeros(nfft,1);      % allocate output spectrum in advance
X_right = zeros(nfft,1);

smoothX_left  = zeros(nfft,1);      % allocate output spectrum in advance
smoothX_right = zeros(nfft,1);

avgAmpX_left  = zeros(nfft,1);
avgAmpX_right = zeros(nfft,1);

LSD = zeros(nframes,1);

SX_left  = zeros(nfft,nframes);      % allocate output spectrogram in advance
SX_right = zeros(nfft,nframes);

% smoothed auto, cross power spectrum and MSC
XC_11 = zeros(nfft,1);        % auto power spectrum
XC_22 = zeros(nfft,1);
XC_12 = zeros(nfft,1);        % cross power spectrum 
MSC = zeros(nfft,1);          % magnitude squared coherence

vad_decision  = zeros(nx,1); 
detect_result = zeros(nx,1);
th_msc = 0.61;    % for pureVoice and voiceMusic
MSC_rate = zeros(nframes,1);

% processed/enhanced audio data
%%%yout = zeros(nframes*nhop+nfft, 2);  % or zeros((nframes-1)*nhop+M, 2); will check???

% for LP based output
yout = zeros(nframes*nhop+M, 2);


times_out = (0:length(yout)-1)*1/Fs;

NN = 20;
kk = 5;

% HP and LP Butterworth filters
Wn = 5000/(Fs/2);                   % 8000 for all before, Normalized cutoff frequency
[Bl,Al] = butter(4,Wn,'low');       % Butterworth filter
%hl = fvtool(Bl,Al);                 % Plot magnitude response
[Bh,Ah] = butter(4,Wn,'high');      % Butterworth filter

% frame processing
for m=0:nframes-1
  xframe_l = x_left_vpu(m*nhop+(1:M));
  xframe_r = x_right_acoustic(m*nhop+(1:M));
  
  % windowing and fft
  xw_l = w .* xframe_l;                         % Apply window
  xw_r = w .* xframe_r;
  
  % calculate LPC coefficients
  A_left  = lpc(xw_l, NN);
  A_right = lpc(xw_r, NN);
  % without windowing
  %A_left  = lpc(xframe_l, NN);
  %A_right = lpc(xframe_r, NN);
  
  x_res = kk * filter(A_left, A_right, xw_l);
  
  % HP for acoustic mic; LP for contact mic
  x_lp = filter(Bl, Al, xw_l);
  x_hp = filter(Bh, Ah, xw_r);
  % Combine LP and HP signals
  x_com = x_lp + x_hp;
  
  %xwzp_l = [xw_l(Mo2+1:M);zp;xw_l(1:Mo2)];      % in zero-phase form!!
  %xwzp_r = [xw_r(Mo2+1:M);zp;xw_r(1:Mo2)];
  xwzp_l = [xw_l;zp];     % to use IFFT directly!!!
  xwzp_r = [xw_r;zp];
  X_left  = fft(xwzp_l);
  X_right = fft(xwzp_r);
  
  % smooth the spectrum
  smoothX_left  = alpha*smoothX_left + (1-alpha)*X_left;
  smoothX_right = alpha*smoothX_right + (1-alpha)*X_right;
  
  % spectrogram of input signal
  SX_left(:,m+1) = X_left;
  SX_right(:,m+1) = X_right;
  
  % coherence  
  XC_11 = alpha*XC_11 + (1-alpha)*(X_left.*conj(X_left));
  XC_22 = alpha*XC_22 + (1-alpha)*(X_right.*conj(X_right));
  XC_12 = alpha*XC_12 + (1-alpha)*(X_left.*conj(X_right));
  % MSC
  MSC = (XC_12.*conj(XC_12)) ./ (XC_11.*XC_22);
  % average of lower frequency 0 - 500Hz, 500/(Fs/nfft) -> 21.xx -> 22
  MSC_ave = sum(MSC(1:44))/44;
  vad_decision((m)*nhop+(1:nhop),1) = 1 - MSC_ave;
  % set flags
  if ((1 - MSC_ave) > th_msc)
     detect_result((m)*nhop+(1:nhop),1) = 0.5;
  else
     detect_result((m)*nhop+(1:nhop),1) = 1.0;
  end
  
  % speech enhancement by equalizer
  
  % average spectrum of before and after j frequency bins, the first and
  % last j values are not averaged
  jj = 5;    %3;
  % without smoothing
%  for mm = jj+1:length(X_left)-jj
%     avgAmpX_left(mm)  = sum(abs(X_left(mm-jj:mm+jj)))/(2*jj+1);
%     avgAmpX_right(mm) = sum(abs(X_right(mm-jj:mm+jj)))/(2*jj+1);
%  end
%  avgAmpX_left  = [abs(X_left(1:jj)); avgAmpX_left(jj+1:length(X_left)-jj); abs(X_left(length(X_left)-jj+1:length(X_left)))];
%  avgAmpX_right = [abs(X_right(1:jj)); avgAmpX_right(jj+1:length(X_right)-jj); abs(X_right(length(X_right)-jj+1:length(X_right)))];
  
  % using smoothed spectrum
  for mm = jj+1:length(smoothX_left)-jj
     avgAmpX_left(mm)  = sum(abs(smoothX_left(mm-jj:mm+jj)))/(2*jj+1);
     avgAmpX_right(mm) = sum(abs(smoothX_right(mm-jj:mm+jj)))/(2*jj+1);
  end
  avgAmpX_left  = [abs(smoothX_left(1:jj)); avgAmpX_left(jj+1:length(smoothX_left)-jj); abs(smoothX_left(length(smoothX_left)-jj+1:length(smoothX_left)))];
  avgAmpX_right = [abs(smoothX_right(1:jj)); avgAmpX_right(jj+1:length(smoothX_right)-jj); abs(smoothX_right(length(smoothX_right)-jj+1:length(smoothX_right)))];
  
  % get amplitude ratio
  Heq = avgAmpX_right ./ avgAmpX_left;
  % get equalized contact mic output
  equalizedX_left = (abs(X_left) .* Heq) .* exp(1j .* angle(X_left));
  
  % calculate Log Spectral Deviation
  LSD(m+1, 1) = sqrt(sum((log10(equalizedX_left.*conj(equalizedX_left))-log10(X_right.*conj(X_right))) .^2));
  
  % processed/enhanced output
%  y_left = real(ifft(X_left));
  %%%y_left = real(ifft(equalizedX_left));
  %%%yout(m*nhop+(1:nfft), 1) = yout(m*nhop+(1:nfft), 1) + y_left;
  %%%y_right = real(ifft(X_right));
  %%%yout(m*nhop+(1:nfft), 2) = yout(m*nhop+(1:nfft), 2) + y_right;
  
  %y_left = real(ifft(equalizedX_left));
  %yout(m*nhop+(1:M), 1) = yout(m*nhop+(1:M), 1) + x_res;
  yout(m*nhop+(1:M), 1) = yout(m*nhop+(1:M), 1) + x_com;
  %y_right = real(ifft(X_right));
  yout(m*nhop+(1:M), 2) = yout(m*nhop+(1:M), 2) + xw_r;

end

  % mean of LSD
  meanLSD = mean(LSD);

  % spectrograms
  t = (M/2:nhop:M/2+(nframes-1)*nhop)/Fs;    % time in seconds
  f = 0.001*(0:nfft/2)*Fs/nfft;              % frequency
  Xdb_l = 20*log10(abs(SX_left));
  Xdb_r = 20*log10(abs(SX_right));
  Xdb_left = Xdb_l(1:nfft/2+1, :);
  Xdb_right = Xdb_r(1:nfft/2+1, :);
    
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
  %title('Time domain voice signal with background music');
  title('Time domain clean voice signals');
  %axis([0 4.0 -1.1 1.1]);
  axis([0 5.0 -1.1 1.1]);
  legend('Acoustic microphone', 'Contact microphone');

  subplot(2,1,2);
  plot(times_out, yout(:,2)', 'b');
  set(gca, 'FontName', 'Arial', 'FontSize', 12);   %'linewidth', 2
  hold on;
  grid on;
  plot(times_out, yout(:,1)', 'y');
  xlabel('Time (Second)');
  ylabel('Amplitude');
  title('Time domain processed/enhanced voice');
  %axis([0 4.0 -1.1 1.1]);
  axis([0 5.0 -1.1 1.1]);
  legend('Acoustic microphone', 'Contact microphone');
  
  figure (2);
  plot(1:nframes, LSD);
  set(gca, 'FontName', 'Arial', 'FontSize', 12);
  hold on;
  grid on;
  plot(1:nframes, meanLSD*ones(nframes,1), 'r', 'linewidth', 2);
  xlabel('Frame Number');
  ylabel('Amplitude');
  title('Log Spectrum Deviation');
  legend('LSD', 'mean of LSD');
  
  % output the processed audiodata
  %audiowrite('enhancedData.wav', yout, Fs);
  %%%audiowrite('enhancedData_clean.wav', yout, Fs);
  
  %audiowrite('enhancedData_HPLP_filter.wav', yout, Fs);
  %audiowrite('enhancedData_HPLP_filter_again.wav', yout, Fs);
  

  %---------------------------------------------
  % 9/4/2021
  %---------------------------------------------
  % OTOG2_case01, set output file name
  %outputFileName = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind3mpsFanB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind3mpsFanB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
  % PI
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind3mpsFanB.wav');
  % PI
  outputFileName_wns = strcat(outputFolderName_ens, 'sc_HpfLpf_eh_HPF.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, inputFile, '\HpfLpf_eh_HPF_wns_OTOG2_SpeechWind3mpsFanB_mic7_1ch.wav');
  % PI
  outputFileName_wns_1ch = strcat(outputFolderName_ens, 'sc_HpfLpf_eh_HPF_1ch.wav');
  
  % OTOG2_case02, set file name 
  %outputFileName = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind3mpsTunnelA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind3mpsTunnelA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind3mpsTunnelA_mic7_1ch.wav');
  
  % OTOG2_case03, set file name 
  %outputFileName = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind5mpsFanA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind5mpsFanA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind5mpsFanA_mic7_1ch.wav');
  
  % OTOG2_case04, set file name 
  %outputFileName = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'HpfLpf_eh_OTOG2_SpeechWind5mpsTunnelB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_wns_OTOG2_SpeechWind5mpsTunnelB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'HpfLpf_eh_HPF_wns_OTOG2_SpeechWind5mpsTunnelB_mic7_1ch.wav');
  
  %---------------------
  % save
  %---------------------
  % without WNS
  %audiowrite(outputFileName, yout, Fs);
  %outInfo = audioinfo(outputFileName)
  % with WNS
  audiowrite(outputFileName_wns, yout, Fs);
  outInfo = audioinfo(outputFileName_wns)
  
  %-----------------------
  % for comparison
  %-----------------------
  % without WNS
  %audiowrite(outputFileName_1ch, yout(:,1), Fs);
  % with WNS
  audiowrite(outputFileName_wns_1ch, yout(:,1), Fs);
  % PI
  disp(outputFolderName_ens);