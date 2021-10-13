  %%%%%%%%%%   LPF-HPF method uding adaptive crossover   %%%%%%%%%%%%%%%%%
  
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

pkg load signal;

clc;
clear all; 
##close all;
addpath('/Users/harry/Desktop/PyLearn/Dev/WindNoiseSupress/ShivaKumarChikkala/ShivaKumarChikkala/windoniseAlgo/Matlab_Code/V3_12Oct')
tic
 %-------------------------------------------------------------------------
  % load input data
  %-------------------------------------------------------------------------
%  if(0)
  [x_in, Fs] = audioread('OTOG2_Speech_silent_wind3_wind5_FanB_mic7.wav');
##  x_in = x_in(1:Fs*10,:);      % take only 8 seconds of input data
%  Fs = 48000;
%  endif
  %-------------------------------------------------------------------------
  % pre treat input data with HPF with cutoff frequency of 80 Hz
  %-------------------------------------------------------------------------

%  x_in = ones(100000,2);
%  p = ones(50000,1); q = zeros(50000,1);
%%  x_in = [p;q];
%  x_in(:,1) = [p;p];
%  x_in(:,2) = [p;q];
%  
%%  x_in = ones(100000,2);
%  x_right_acoustic = x_in(:,1);
%  x_left_vpu = x_in(:,2);
 
  [b, a]            = butter(6, 80/(48e3/2), 'high');     % generate filter coeffcients using butterworth filter of 6 order
  x_right_acoustic  = filter(b, a, x_in(:,1));            % filter acoustic mic data with above coeffcients
  x_left_vpu        = filter(b, a, x_in(:,2));            % filter contact mic data with above coeffcients

  %-------------------------------------------------------------------------
  % Filter deisgn for band-pass filters and pre-treat time data for contact mic only (to estimate noise thresholds)
  %-------------------------------------------------------------------------
 
  band_freqs=([100 1500]);
  
  % band-pass filter      % investigate order of filter
  
%  [b_high,a_high]     = butter(2,band_freqs(1)/Fs,'high');
%  [b_low,a_low]       = butter(2,band_freqs(2)/Fs,'low');
%  b2                  = conv(b_high,b_low);
%  a2                  = conv(a_high,a_low);
%  x_left_vpu_filtered = filter(b2,a2,x_in(:,2));
  
  [b_band,a_band]     = butter(2,[band_freqs(1)/Fs, band_freqs(2)/Fs]);
  x_left_vpu_filtered = filter(b_band, a_band, x_in(:,2));
  
  %--------------------------------------------------------------------------
  % start algorithm
  %--------------------------------------------------------------------------
  
%  times           = (0:length(x_left_vpu)-1)*1/Fs;          % get time range values in seconds  % just for time domain plot
  nx = length(x_left_vpu);      % assume length(x) > M
  noverlap        = 512;                                    % set up signal analysis parameters
  window          = hann(1024, 'periodic');                      
  M               = length(window);
  w               = window(:);   
  nhop            = M-noverlap;
  nx              = length(x_left_vpu);                     % assume length(x) > M
  nframes         = 1+fix((nx-M)/nhop);                     % another way to calculate frame numbers
  xframe_l        = zeros(M,1);                             % framed input samples
  xframe_r        = zeros(M,1);
  nfft=2048;
  Modd = mod(M,2);   % 0 if M even, 1 if odd
  Mo2 = (M-Modd)/2;  % Make sure it's a column 
  zp = zeros(nfft-M,1);         % zero-padding for each FFT
  XC_11 = zeros(nfft,1);        % auto power spectrum
  XC_22 = zeros(nfft,1);
  XC_12 = zeros(nfft,1);        % cross power spectrum 
  MSC = zeros(nfft,1); 
  % RMS diff
  rms_diff = zeros(M,1);
  RMSD_decision = zeros(nx,1);
  RMSD_rate = zeros(nframes,1);
  
  exp_avg_r       = 0;
  NE              = zeros(M,1);
%  exp_avg_r_out   = zeros(nframes*nhop+M, 1);
  exp_avg_r_out   = zeros(nframes, 1);
%  sum_vad_avg     = zeros(nframes*nhop+M, 1);
  times           = (0:length(nframes)-1)*1/Fs;
  yout            = zeros(nframes*nhop+M, 2);
  times_out       = (0:length(yout)-1)*1/Fs;
  filtered_output = zeros(nframes*nhop+M, 2);
  VAD_signal      = zeros(M,1);
  
  VAD_delay1 = 0; 
  VAD_delay = 0; 
  E_out           = zeros(M,1);
  vad_t_update           = zeros(M,1); 
  VAD_signal1     = zeros(nframes*nhop+M, 1);
  SpeechFlag     = zeros(nframes*nhop+M, 1);
%  vad_threshold   = 2e-3;
  vad_threshold   = 2e-5;
  counter =18
  E_full          =  zeros(nframes, 1);         % edit this to proper length
  exp_avg_r_out_db =  zeros(nframes*nhop+M, 1);
  sum_ne          = zeros(M,1);
  Wn              = ones(nframes,1);
  
  %--------------------------------------------------------------------------
  %  Cutoff frequencies for HP and LP Butterworth filters
  %--------------------------------------------------------------------------
  
  Xc        = [ 80, 315, 630, 1250, 2500, 5000];           % frequencies to define crossover

  %--------------------------------------------------------------------------
  %         frame based processing
  %--------------------------------------------------------------------------

  %% VAD Parameters
%  tA = 1/(Fs *0.01); % 0.0025 s   % Time constants of the VAD
  tA = 0.01;  % Time constants of the VAD
%  tA = 0;  % Time constants of the VAD
  % SNR detection threshold
%  SNR_THS = 1;
%  K = 256;
  E = 0;
%  Eflr = 0;
%   x_in = ones(100000,2);

  alpha = 0.91;  % smooth factor 0.8 for Lin method
  E_for_all_Samples = zeros(nx,1);
  vadthr_all_Samples = zeros(nx,1);

  for m = 0:nframes - 1
 
   xframe_l           = x_left_vpu(m*nhop+(1:M));
   xframe_l_filtered  = x_left_vpu_filtered(m*nhop+(1:M));
   xframe_r           = x_right_acoustic(m*nhop+(1:M));

   
  %--------------------------------------------------------------------------
  %                     VAD
  %--------------------------------------------------------------------------

%    counter =18
%    if(VAD_delay > 0.7) % can change freom 0.7,0.8,1
%    counter = 0;
%    end
%
%    if(counter < 18)
%    VAD_signal = ones(M,1);
%    counter = counter +1;
%    continue
%  else
%    counter = 0;
%    [VAD_signal, E_out, VAD_delay1] = VAD(xframe_l_filtered, E, tA, M, nhop, vad_threshold, Fs, VAD_delay);
%        E = E_out(nhop);      % update E with previous average out value
%        VAD_delay = VAD_delay1;
%        VAD_signal1(m*nhop+(1:nhop), 1) = VAD_signal(1:nhop);
%        E_full(m*nhop+(1:nhop),1)       = E_out(nhop);
%     end
  
%  if(0)
alphaforVad_t =0.1;
    
    [VAD_signal, E_out] = VAD(xframe_l_filtered, E, tA, M, nhop, vad_threshold, Fs);
    if(E_out(nhop) < E)
        vad_threshold = vad_threshold + alphaforVad_t*E;
     else
        vad_threshold = 2e-5;
     endif
      if(  vad_threshold <= 2e-5)
          vad_threshold   = 2e-5;
        endif
    E = E_out(nhop);      % update E with previous average out value

    E_for_all_Samples(m*nhop+(1:M)) = E + zeros(M,1);
    vadthr_all_Samples(m*nhop+(1:M)) = vad_threshold + zeros(M,1);
%    VAD_delay = VAD_delay1;
    VAD_signal1(m*nhop+(1:nhop), 1) = VAD_signal(1:nhop);
    E_full(m*nhop+(1:nhop),1)       = E_out(nhop);
%    endif
  %--------------------------------------------------------------------------
  %                     VAD by Lin
  %--------------------------------------------------------------------------
  if(0)
%   [detect_result, AllT_rate]  = vad_rmsd(x_left_vpu, x_right_acoustic, Fs, m); 
  
  % VAD_02: SSC
  xw_l = w .* xframe_l;                         % Apply window
  xw_r = w .* xframe_r;
  xwzp_l = [xw_l(Mo2+1:M);zp;xw_l(1:Mo2)];      % in zero-phase form!!
  xwzp_r = [xw_r(Mo2+1:M);zp;xw_r(1:Mo2)];
  %%%xwzp = [xw;zp];
  %xwzp_l = [xw_l;zp];
  %xwzp_r = [xw_r;zp];
  X_left  = fft(xwzp_l);
  X_right = fft(xwzp_r);
  
  % spectrogram of input signal
  SX_left(:,m+1) = X_left;
  SX_right(:,m+1) = X_right;
  
  XC_11 = alpha*XC_11 + (1-alpha)*(X_left.*conj(X_left));
  XC_22 = alpha*XC_22 + (1-alpha)*(X_right.*conj(X_right));
  XC_12 = alpha*XC_12 + (1-alpha)*(X_left.*conj(X_right));
  SSC_n_1 = 0;
  SSC_d_1 = 0;
  SSC_n_2 = 0;
  SSC_d_2 = 0;
  for mm = 1:40          % 20 for windTest signal
      SSC_n_1 = SSC_n_1 + XC_11(mm)*mm*Fs/nfft;
      SSC_d_1 = SSC_d_1 + XC_11(mm);
      SSC_n_2 = SSC_n_2 + XC_22(mm)*mm*Fs/nfft;
      SSC_d_2 = SSC_d_2 + XC_22(mm);
  end
  SSC_1      = SSC_n_1 / SSC_d_1;
  SSC_1_flag = ((40*Fs/nfft) - SSC_1) / (40*Fs/nfft);     % 20 for windTest signal
  SSC_decision_1((m)*nhop+(1:nhop),1) = SSC_1_flag;
  SSC_2      = SSC_n_2 / SSC_d_2;
  SSC_2_flag = ((40*Fs/nfft) - SSC_2) / (40*Fs/nfft);
  SSC_decision_2((m)*nhop+(1:nhop),1) = SSC_2_flag;
  
  SSC_SX(:,m+1) = ones(nfft,1)*SSC_1_flag;

  % use VPU only
  %if (SSC_1_flag > 0.7)     % 0.7 for windVoice and pureVoice
  if (SSC_1_flag > 0.65)    % 0.65 for voiceMusic   
      SSC_rate(m+1,1) = 0;
%      detect_result((m)*nhop+(1:nhop),1) = 0.5;    
  else
      SSC_rate(m+1,1) = 1;
%      detect_result((m)*nhop+(1:nhop),1) = 1.0;
  end
  
  %% RMS_DIFF
  
  rms_left    = rms(xframe_l);
  rms_dB_left = 20*log10(rms_left+1.0e-12);
  rms_right    = rms(xframe_r);
  rms_dB_right = 20*log10(rms_right+1.0e-12);
  
  rms_diff = 0.95*rms_diff + (1-0.95)*(rms_dB_right - rms_dB_left);  % 0.98 for windVoice
  RMSD_decision((m)*nhop+(1:nhop),1) = mean(rms_diff);
  
  RMSD_SX(:,m+1) = ones(nfft,1)*mean(rms_diff);

  % use RMS diff
  %if (mean(rms_diff) > 18.5)      % 18.5 for windVoice signal
%  if (mean(rms_diff) > -1.0)      % -1.0 for pureVoice signal
  if (mean(rms_diff) > 18.25)      % 18 for VoiceMusic signal    
      RMSD_rate(m+1,1) = 0;
%      detect_result((m)*nhop+(1:nhop),1) = 0.5;    
  else
      RMSD_rate(m+1,1) = 1;
%      detect_result((m)*nhop+(1:nhop),1) = 1.0;
  end

  % combined RMSD and SSC_VPU
  if (SSC_1_flag > 0.65) || (mean(rms_diff) > 18.25)   % 0.65 for voiceMusic    
%      detect_result((m)*nhop+(1:nhop),1) = 0.5;
      detect_result((m)*nhop+(1:nhop),1) = 0;
      AllT_rate(m+1,1) = 0;
  else
      detect_result((m)*nhop+(1:nhop),1) = 1.0;
      AllT_rate(m+1,1) = 1;
  end

  %--------------------------------------------------------------------------
  %           Adaptive crossover frequencies
  %--------------------------------------------------------------------------
  
##   [x_com, xw_l, xw_r, NE, sum_ne, Wn]  =  adaptive_crossover(xframe_l, xframe_r, Xc, Fs, exp_avg_r, VAD_signal, M, w, sum_ne);
##    exp_avg_r = NE(nhop);
##    exp_avg_r_out(m*nhop+(1:nhop), 1)     = NE(1:nhop);
##%    sum_vad_avg += NE(1:nhop);
##    exp_avg_r_out_db = 10*log10(exp_avg_r_out+1.0e-12);
##%     sum_vad_avg(m*nhop+(1:M), 1) = sum_ne;
##%    sum_vad_avg(m*nhop+(1:M), 1) = disp(sum(NE));
##    yout(m*nhop+(1:M), 1)   = yout(m*nhop+(1:M), 1) + x_com;
##    yout(m*nhop+(1:M), 2)   = yout(m*nhop+(1:M), 2) + xw_r;
##    Wn(m+1,1) = Wn;
%   disp(Wn);
  endif % end of if(0)
  end
       
%**************************************************************************
%       Present and plot the VAD decision results overlayed the test sequence.
%**************************************************************************
%if(0)
%    t = (0:length(s)-1)/fs;
    figure;
%    subplot(2,1,1)
    plot(x_left_vpu./674.35, 'b');
    hold on; plot( (VAD_signal1./674.35)', 'r'); grid;
    hold on; plot(E_for_all_Samples, 'k');
    hold on; plot(vadthr_all_Samples, 'g');
%        hold on; plot( VAD_signal1', 'r', 'linewidth', 1.5);
    title('VAD with time constant 0.01 with noisy speech'); 
    legend('contact mic', 'VAD decision (1=speech, 0=noise)','Energy','vad_threshold');
%    hold on; plot(exp_avg_r_out_db,'k'); 
%    subplot(2,1,2)
%    plot(x_left_vpu, 'b'); hold on;
%    plot(detect_result, 'r');
%    title(' Lins VAD with noisy speech'); 
%    xlabel('samples');  ylabel('Signal value');
%    legend('contact mic', 'Lins VAD decision (1=speech, 0=noise)');
%    grid on;

%  figure; plot(x_left_vpu, 'b'); hold on; plot( VAD_signal1', 'r'); grid;
  
%    axis([t(1), t(end), -0.5, 1.5]);
%end
    
%    figure; plot(x_left_vpu); hold on; plot(10*log10(E_full+1.0e-12)); grid; legend('input', 'vad avg'); title('exp vad avg with time constant of 0.025 seconds ');
%    figure; plot(x_right_acoustic); hold on; plot(exp_avg_r_out_db,'k'); hold on; grid; legend('input', 'noise estimate'); title('noise est-Xover with time constant of 0,01 seconds'); 
    
%    elapsed_time = toc  % stop timer

%    audiowrite('C:\Users\00317859\OneDrive - Sigma AB\Skrivbordet\OTO\copy of 2021-09-15 From Lin - sandbox\outputFiles\LPF-HPF_crossover_output.wav', yout, Fs);
%    audiowrite('C:\Users\00317859\OneDrive - Sigma AB\Skrivbordet\OTO\copy of 2021-09-15 From Lin - sandbox\outputFiles\LPF-HPF_crossover_output_1ch.wav', yout(:,1), Fs);
%  end  
%  end
toc