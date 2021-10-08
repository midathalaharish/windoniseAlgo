%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Speech enhancement by equalizer or FIR filter.
% Reference paper: Reconstruction Filter Design for Bone-Conducted Speech    
%                  ON EQUALIZATION OF BONE CONDUCTED SPEECH FOR IMPROVEDSPEECH QUALITY
%                  Distortion Measures for Speech Processing
% Input:       x_left_vpu (channel-1): contact mic signal  
%              x_right_acoustic (channel-2): acoustic mic signal
%
% Output:      yout - two channel data file
%                     channel-1 (left):  enhanced vpu signal
%                     channel-2 (right): ref acousitc signal
%
% Test files:  iTP02 - rec_alexa_pizza.wav and 
%              iTP02_Music_rec_alexa_pizza.wav are from Sigma
%
% Date:        4/19/2021, initial, based on VAD_Algorithm_final.m
% Modified:    4/21/2021, add smooth filter for spectrum and distortion 
%                         measure: Log Spectral Deviation.
%
% Notes:       5/05/2021, need to do more test and comparisons
% Modified:    8/03/2021, tested with OTO gen2 data.
% Modified:    8/14/2021, tested with OTO gen2 data with 4 directions of wind
% Modified:    8/28/2021, verify OTOG2 new data from Gongqiang, add common
%                         folderName and fileName. Create 1-ch data for
%                         comparison
%                         Applied WNS to input file, then run SE
% Modified:    9/03/2021, continue on OTOG2 data
% Modified:    9/04/2021, add HPF for the acoustic mic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean up
clc;
clear; 
close all;

%-------------------------------------------------------------------------
% load input data
%-------------------------------------------------------------------------

%folderName = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\';
% PI
folderName = 'C:\Users\00317859\Desktop\OTO\copy of 2021-09-15 From Lin - sandbox\inputFiles\';

%%outputFolderName = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\EnhancedData\';  % SC
% PI
outputFolderName = 'C:\Users\00317859\Desktop\OTO\copy of 2021-09-15 From Lin - sandbox\outputFiles\';

%folderName_wns = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\WindSuppressed\';
% PI
folderName_wns = 'C:\Users\00317859\Desktop\OTO\copy of 2021-09-15 From Lin - sandbox\inputFiles\';

%outputFolderName_ens = 'C:\Projects\AfterAug2\OTOG2_Data_fromGongqiang\CreatedData\WindSuppressed\EnhancedData\';
% PI
outputFolderName_ens = 'C:\Users\00317859\Desktop\OTO\copy of 2021-09-15 From Lin - sandbox\outputFiles\';


% load recorded data from Sigma, a two-channel wave file sampling @ 48 KHz, 
% first, left channel is contact mic signal, second, right channel is 
% acoustic mic

%[x_in, Fs] = audioread('iTP02_Music_rec_alexa_pizza.wav');    % speech+music
%[x_in, Fs] = audioread('iTP02_rec_alexa_pizza.wav');          % clean speech

% get contact mic VPU signal and acoustic mic signal
%x_left_vpu = x_in(:,1);
%x_right_acoustic = x_in(:,2);

% load recorded data from Sage, a 18-channel wave file sampling @ 48 KHz.
% channel 1, Mouth Mic, channel 2, Stethoscope, 
% channel -3-18, OTO mics, where channels 9-10 contact mics.

% Case0110: speech_noise harvard_sentence normal wind 3m/s Front
%           We talked of the slide show in the circus.
%[x_in, Fs] = audioread('case0110_session05_09_2021_09_09_54_dv1.wav');

% Case0158: speech_noise harvard_sentence normal wind 5m/s Front
% A tame squirrel makes a nice pet.
%[x_in, Fs] = audioread('case0158_session05_09_2021_09_09_54_dv1.wav');

% Case0208: speech_noise harvard_sentence normal wind 8m/s Front
%           hey facebook whats a recipe for pizza.
%[x_in, Fs] = audioread('case0208_session05_09_2021_09_09_54_dv1.wav');

% get contact mic VPU signal and acoustic mic signal
%x_left_vpu = x_in(:,9)*50;           %mic 7    as original data are weak so x50 about 33 dB gain
%x_right_acoustic = x_in(:,17)*50;    %mic 15

%--------------------------------------------------
% 8/14/2021 for 4 directions of wind
%--------------------------------------------------
%-----------------------------------------------------------------------
% Case 01:  For wake word normal + wind 8m/s Front
%           hey facebook whats a recipe for pizza.
% Choose Mic16 (better than Mic15 by listening) and VPU7
%----------------------------------------------------------------------- 
% load mic16 data as input x_right_acoustic
%[x_right_acoustic, Fs] = audioread('FrontWind8mps_acoustic16.wav');
%x_right_acoustic = x_right_acoustic * 25;  % as original data are weak so x50 about 33 dB gain
% load mic7 vpu data as input x_left_vpu   % x50 get clipping, so change to x25. 
%[x_left_vpu, Fsv] = audioread('FrontWind8mps_vpu7.wav');
%x_left_vpu = x_left_vpu * 25;

%-----------------------------------------------------------------------
% Case 02:  For wake word normal + wind 8m/s Left
%           Hey Facebook, what's the weather in Cleveland today?
% Choose Mic16 and VPU7
%----------------------------------------------------------------------- 
% load mic16 data as input x_right_acoustic
%[x_right_acoustic, Fs] = audioread('LeftWind8mps_acoustic16.wav');
%x_right_acoustic = x_right_acoustic * 25;  % as original data are weak so x50 about 33 dB gain
% load mic7 vpu data as input x_left_vpu   % x50 get clipping, so change to x25. 
%[x_left_vpu, Fsv] = audioread('LeftWind8mps_vpu7.wav');
%x_left_vpu = x_left_vpu * 25;

%-----------------------------------------------------------------------
% Case 03:  For wake word normal + wind 8m/s Back
%           Hey Facebook, what sound does a turtle make?
% Choose Mic15 (similar to Mic16) and VPU7
%-----------------------------------------------------------------------
% load mic15 data as input x_right_acoustic
%[x_right_acoustic, Fs] = audioread('BackWind8mps_acoustic15.wav');
%x_right_acoustic = x_right_acoustic * 25;
% load mic7 vpu data as input x_left_vpu
%[x_left_vpu, Fsv] = audioread('BackWind8mps_vpu7.wav');
%x_left_vpu = x_left_vpu * 25;

%-----------------------------------------------------------------------
% Case 04:  For wake word normal + wind 8m/s Right
%           Hey Facebook, send a message to Claire
% Choose Mic15 and VPU7
%-----------------------------------------------------------------------
% load mic15 data as x_right_acoustic
%[x_right_acoustic, Fs] = audioread('RightWind8mps_acoustic15.wav');
%inputInfo1 = audioinfo('RightWind8mps_acoustic15.wav')
%x_right_acoustic = x_right_acoustic * 25;
% load mic7 vpu data as x_left_vpu
%[x_left_vpu, Fsv] = audioread('RightWind8mps_vpu7.wav');
%inputInfo2 = audioinfo('RightWind8mps_vpu7.wav')
%x_left_vpu = x_left_vpu * 25;

%--------------------------------------------------
% 8/28/2021 for new OTOG2 data
%--------------------------------------------------
%-----------------------------------------------------------------------
% OTOG2_case01: speech_wind_3mps_FanB, the input file is created by
%               CreateTwoChannelData.m
%               ch-1: acoustic mic; ch-2: contact mic
%-----------------------------------------------------------------------
% set file name
%fileName = strcat(folderName, 'OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
%fileName_wns = strcat(folderName_wns, 'wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');

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
% PI
fileName_wns = strcat(folderName_wns, 'A007_OTOG2_SpeechWind3mpsFanB_mic15_speechwind3mps_mic7.wav');

% get data from the file without WNS
%[x_in, Fs]  = audioread(fileName);
%inputInfo   = audioinfo(fileName)

% get data from the file with WNS
[x_in, Fs]  = audioread(fileName_wns);
audioinfo(fileName_wns)

% add HPF for acoustic mic
[b, a] = butter(6, 300/(48e3/2), 'high');     %9/4/2021

x_left_vpu = x_in(:,2);
%x_right_acoustic = x_in(:,1);
x_right_acoustic = filter(b, a, x_in(:,1));   %9/4/2021

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
yout = zeros(nframes*nhop+nfft, 2);  % or zeros((nframes-1)*nhop+M, 2); will check???

times_out = (0:length(yout)-1)*1/Fs;

% frame processing
for m=0:nframes-1
  xframe_l = x_left_vpu(m*nhop+(1:M));
  xframe_r = x_right_acoustic(m*nhop+(1:M));
  
  % windowing and fft
  xw_l = w .* xframe_l;                         % Apply window
  xw_r = w .* xframe_r;
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
  %% The magnitude-squared coherence estimate is a function of frequency with values between 0 and 1. 
  %%These values indicate how well x corresponds to y at each frequency. The magnitude-squared coherence is a function of the 
  %%power spectral densities, Pxx(f) and Pyy(f), and the cross power spectral density, Pxy(f), of x and y:
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
  y_left = real(ifft(equalizedX_left));
  yout(m*nhop+(1:nfft), 1) = yout(m*nhop+(1:nfft), 1) + y_left;
  y_right = real(ifft(X_right));
  yout(m*nhop+(1:nfft), 2) = yout(m*nhop+(1:nfft), 2) + y_right;

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
    
  %------------------------------------------------------------------------
  % set up figure
  %------------------------------------------------------------------------
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
  
  %------------------------------------------------------------------------
  % output the processed audiodata:
  % channel-1 (left):  enhanced vpu signal
  % channel-2 (right): ref acousitc signal
  %------------------------------------------------------------------------
  
  %audiowrite('enhancedData.wav', yout, Fs);
  %audiowrite('enhancedData_clean.wav', yout, Fs);
  %audiowrite('enhancedData_music.wav', yout, Fs);
  
  %audiowrite('eh_normalSpeech_wind3mps.wav', yout, Fs);
  %audiowrite('eh_normalSpeech_wind5mps.wav', yout, Fs);
  %audiowrite('eh_normalSpeech_wind8mps.wav', yout, Fs);
  
  %audiowrite('eh_FrontWind8mps.wav', yout, Fs);
  %audiowrite('eh_LeftWind8mps.wav', yout, Fs);
  %audiowrite('eh_BackWind8mps.wav', yout, Fs);
  %audiowrite('eh_RightWind8mps.wav', yout, Fs);

  % OTOG2_case01, set output file name
  %outputFileName = strcat(outputFolderName, 'eh_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'eh_OTOG2_SpeechWind3mpsFanB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind3mpsFanB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind3mpsFanB_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind3mpsFanB_mic7_1ch.wav');
  
  % OTOG2_case02, set file name 
  %outputFileName = strcat(outputFolderName, 'eh_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'eh_OTOG2_SpeechWind3mpsTunnelA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind3mpsTunnelA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind3mpsTunnelA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind3mpsTunnelA_mic7_1ch.wav');
  
  % OTOG2_case03, set file name 
  %outputFileName = strcat(outputFolderName, 'eh_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'eh_OTOG2_SpeechWind5mpsFanA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind5mpsFanA_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind5mpsFanA_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind5mpsFanA_mic7_1ch.wav');
  
  % OTOG2_case04, set file name 
  %outputFileName = strcat(outputFolderName, 'eh_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
  %outputFileName_1ch = strcat(outputFolderName, 'eh_OTOG2_SpeechWind5mpsTunnelB_mic7_1ch.wav');
  %outputFileName_wns = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
  %outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_wns_OTOG2_SpeechWind5mpsTunnelB_mic7_1ch.wav');
  outputFileName_wns = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind5mpsTunnelB_mic15_mic7.wav');
  outputFileName_wns_1ch = strcat(outputFolderName_ens, 'eh_HPF_wns_OTOG2_SpeechWind5mpsTunnelB_mic7_1ch.wav');
  
  
  % save
  %audiowrite(outputFileName, yout, Fs);
  %outInfo = audioinfo(outputFileName)
  audiowrite(outputFileName_wns, yout, Fs);
  outInfo = audioinfo(outputFileName_wns)
  
  % for comparison
  %audiowrite(outputFileName_1ch, yout(:,1), Fs);
  audiowrite(outputFileName_wns_1ch, yout(:,1), Fs);
  