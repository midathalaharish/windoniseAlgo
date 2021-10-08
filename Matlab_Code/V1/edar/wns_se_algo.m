function [y, R, Sout] = wns_se_algo(x, fs, S)
% [y, R, Sout] = wns_se_algo(x, fs, S)
% S = wns_se_algo   (no input arguments)
%
% (Function description below not updated!!!)
% INPUT SIGNALS
%
% x               sensor signals in columns    [x_vpu x_mic]
%
% INPUT PARAMETERS
%
% fs              sample rate [Hz]
% nfft            number of points in FFT, NFFT.
% noverlap        number of samples the sections of A overlap. If negative,
%                 -NOVERLAP is the "hop size". (The overlap is the window
%                 length minus the hop size.) NOVERLAP < M
% nwin            length of window
% alpha           time weighting factor for spectrum and coherence
% fhp             cutoff frequency of HP filter
% th_msc          threshold value for VAD
%
% OUTPUT SIGNAL
%
% y               processed signal and reference signal in two-column array
%
% OUTPUT DATA
%
% SX_vpu          STFT's of input signal
% SX_mic          STFT's of input signal
% vad_decision
% detect_result
% LSD             Log Spectral Deviation


%--------------------------------------------------------------------------
% Signal analysis parameters
%--------------------------------------------------------------------------
% Default parameter values
if nargin<3
    % Frames processing parameters
    S.nwin                 = 512 ;   % window (frame) length
    S.noverlap             = 256 ;   % frame (window) overlap
    S.nfft                 = 2048 ;
    
    % High-pass filter (acoustic mic) parameters
    S.hp_f                 = 300 ;   % [Hz]
    
    % Wind noise suppression general parameters
    S.WNS_method           = 'RMSD' ; % 'MSC'|'SSC'|'RMSD'|'none'
    S.WNS_counter_th       = 10 ;    % WNS "attack time" [frames]
    S.WNS_frameAttenuation = 0.15;   %0.1
    
    % MSC (magnitude squared coherence) parameters
    S.MSC_alpha            = 0.91;   % psd time weighting factor % 0.8
    S.MSC_f                = 500 ;   % [Hz]
    S.MSC_value_th         = 0.61;   % MSC threshold
    
    % SSC parameters
    S.SSC_f                = 450 ;   % [Hz]
    S.th_ssc               = 0.55;   % SSC threshold % 0.7, 0.4
    
    % RMSD parameters
    S.RMSD_alpha           = 0.95 ;  % 0.98, 0.85 for Long file
    
    % Speech enhancement general parameters
    S.SE_method            = 'HPLP' ; % 'STEQ'|'LP'|'HPLP'|'none'
    
    % STEQ
    S.STEQ_alpha_X         = 0.1;    %0.5; STFT time weighting factor
    S.STEQ_alpha_XC        = 0.1;    %0.5; psd time weighting factor
    
    % LP
    S.LP_order             = 20 ;    % order of linear predictor
    
    % HPLP
    S.HPLP_fx              = 5000 ;  % cross-over frequency [Hz]
    
end

%--------------------------------------------------------------------------
% Trick to provide the default settings struct S
if nargin==0
    y = S ;
    return
end



%--------------------------------------------------------------------------
% Frame processing parameters
S.nwin   = S.nwin - rem(S.nwin,2) ; % nwin must be even (for practical reasons)
nwin     = S.nwin ;
S.noverlap = min(S.noverlap,nwin-1) ; % noverlap must be smaller than nwin
noverlap = S.noverlap ;
S.nfft   = max(S.nfft,nwin) ; % nfft cannot be smaller than nwin
nfft     = S.nfft ;
win      = hann(S.nwin, 'periodic') ;
nhop     = nwin-noverlap ;
nx       = size(x,1) ;      % assume length(x) > nwin
nframes  = floor((nx-nwin)/nhop) ;
% % another way to calculate frame numbers
% nframes = 1+fix((nx-nwin)/nhop);
zp       = zeros(nfft-nwin,1) ;         % zero-padding for each FFT

%--------------------------------------------------------------------------
% MSC (magnitude squared coherence) parameters
i_MSC_f = interp1(  (0:(nfft/2-1))*fs/nfft, 1:nfft/2, S.MSC_f, 'nearest') ; % probably hard-coded in c
S.MSC_f = interp1(  1:nfft/2, (0:(nfft/2-1))*fs/nfft, i_MSC_f, 'nearest') ; % change to actual bin frequency

%--------------------------------------------------------------------------
% SSC parameters
i_SSC_f = interp1(  (0:(nfft/2-1))*fs/nfft, 1:nfft/2, S.SSC_f, 'nearest') ; % probably hard-coded in c
S.SSC_f = interp1(  1:nfft/2, (0:(nfft/2-1))*fs/nfft, i_SSC_f, 'nearest') ; % change to actual bin frequency

%--------------------------------------------------------------------------
% RMSD parameters

%--------------------------------------------------------------------------
% STEQ parameters

%--------------------------------------------------------------------------
% LP parameters (linear prediction filter)

%--------------------------------------------------------------------------
% HPLP parameters (sensor cross-over frequency)
Wn = S.HPLP_fx/(fs/2) ;             % Normalized crossover frequency
[Bl,Al] = butter(4,Wn,'low');       % Butterworth filter
[Bh,Ah] = butter(4,Wn,'high');      % Butterworth filter



%--------------------------------------------------------------------------
% Preallocation of variables (not all needed by MATLAB)
%--------------------------------------------------------------------------
X_vpu = zeros(nfft,1);             % vpu STFT
X_mic = zeros(nfft,1);             % mic STFT

Xtw_vpu = zeros(nfft,1);           % vpu time weighted STFT
Xtw_mic = zeros(nfft,1);           % mic time weighted STFT

avgAmpX_vpu = zeros(nfft,1);       % initial value for spectral smoothing
avgAmpX_mic = zeros(nfft,1);       % initial value for spectral smoothing

R.LSD = nan(nframes,1);

% SX_vpu = zeros(nfft,nframes);      % spectrogram (only for saving STFT's -
% % not needed for real-time processing)
% SX_mic = zeros(nfft,nframes);

% smoothed auto, cross power spectrum and MSC
XC_11 = zeros(nfft,1);             % vpu power spectrum, initial value for time weighting
XC_22 = zeros(nfft,1);             % mic power spectrum, initial value for time weighting
XC_12 = zeros(nfft,1);             % cross power spectrum, initial value for time weighting
%MSC   = zeros(nfft,1);             % magnitude squared coherence

R.MSC_value = nan(nframes,1) ;
R.MSC_wind  = nan(nframes,1) ;

rms_diff = 0 ;                     % initial value for time weighting

% processed/enhanced audio data
y = zeros(nframes*nhop+nfft, 2);   % or zeros((nframes-1)*nhop+M, 2); will check???



%--------------------------------------------------------------------------
% Add HPF for acoustic mic (wind noise rejecting first stage in signal chain)
% *** Make this frame based also!
% *** FB processing order:  WNS -> HP -> SE
%--------------------------------------------------------------------------
[HPmic_b, HPmic_a] = butter(6, S.hp_f/(fs/2), 'high');
%x_mic = filter(b, a, x_mic);
micHPfilterInFrame = false ;
if micHPfilterInFrame
    x(:,2) = filter(HPmic_b, HPmic_a, x(:,2));
end


%--------------------------------------------------------------------------
% Frame processing
%--------------------------------------------------------------------------
WNS_counter = 0 ;

fProgress = waitbar(0,['Estimated remaining time: -' ' s.']) ;
tic
for iframe = 0:nframes-1
    
%     if rem(iframe,1000)==0
%         disp(['iframe = ' num2str(iframe) '/' num2str(nframes)])
%     end
    
    Ihop   = iframe*nhop+(1:nhop) ; % x index range for hop (frame increment)
    Iframe = iframe*nhop+(1:nwin) ; % x index range for (complete) frame
    
    %----------------------------------------------------------------------
    % Get frame data
    %----------------------------------------------------------------------
    xframe_vpu = x(Iframe,1);
    xframe_mic = x(Iframe,2);
    
    %----------------------------------------------------------------------
    % STFT
    %----------------------------------------------------------------------
    % Apply window and fft
    xw_vpu = win .* xframe_vpu;                         % Apply window
    xw_mic = win .* xframe_mic;
    Mo2 = floor(nwin/2) ;
    xwzp_vpu = [xw_vpu(Mo2+1:nwin);zp;xw_vpu(1:Mo2)];      % in zero-phase form!!
    xwzp_mic = [xw_mic(Mo2+1:nwin);zp;xw_mic(1:Mo2)];
    %xwzp_vpu = [xw_vpu;zp];     % to use IFFT directly!!!
    %xwzp_mic = [xw_mic;zp];
    X_vpu = fft(xwzp_vpu);
    X_mic = fft(xwzp_mic);
    
    %     % STFT of input signal (NOTE: No time weighting)
    %     SX_vpu(:,iframe+1) = X_vpu;
    %     SX_mic(:,iframe+1) = X_mic;
    
    %----------------------------------------------------------------------
    % MSC WNS method
    %
    % Parameter  alpha_XC, f_msc, th_msc
    % State      XC_11, XC_22, XC_12 (time weighted psd's and cross-spectrum)
    % Input      X_vpu, X_mic
    % Results    MSC_, wind_decision
    % Output     flag
    %
    % CH: Better to use two acoustic microphones?
    % CH: Better to use equalized sensor signals? Can be a priori weighting or
    % adaptive (estimate when speech in quite conditions is detected)
    %----------------------------------------------------------------------
    % PSD's and cross spectrum (with exponential time weighting)
    XC_11 = S.MSC_alpha*XC_11 + (1-S.MSC_alpha)*(X_vpu.*conj(X_vpu));
    XC_22 = S.MSC_alpha*XC_22 + (1-S.MSC_alpha)*(X_mic.*conj(X_mic));
    XC_12 = S.MSC_alpha*XC_12 + (1-S.MSC_alpha)*(X_vpu.*conj(X_mic));
%     alpha = S.MSC_alpha ;
%     XC_11(1:i_MSC_f) = alpha*XC_11(1:i_MSC_f) + (1-alpha)*(X_vpu(1:i_MSC_f).*conj(X_vpu(1:i_MSC_f)));
%     XC_22 = alpha*XC_22 + (1-alpha)*(X_mic.*conj(X_mic));
%     XC_12 = alpha*XC_12 + (1-alpha)*(X_vpu.*conj(X_mic));
    
    % MSC (magnitude squared coherence). Large MSC => signal is likely not wind
    MSC = (XC_12.*conj(XC_12)) ./ (XC_11.*XC_22);
    
    R.MSC_SX(:,iframe+1) = MSC ; % Result (not to be ported)
    
    % average of lower frequency range, 0 to f_msc
    % CH: (TBI) Use only frequencies where both sensors have sufficient SNR in quiet conditions
    %MSC_ave = mean(MSC(1:i_MSC_f)) ; % scalar, i_MSC_f is index of element at f_msc
    MSC_ave = sum(MSC(1:i_MSC_f))/i_MSC_f ; % scalar, i_MSC_f is index of element at f_msc
    MSC_value = 1 - MSC_ave ;
    % classification
    if (MSC_value > S.MSC_value_th)
        MSC_wind = 1 ;  % Classify as wind
    else
        MSC_wind = 0 ;  % (not wind)
    end
    
    R.MSC_value(iframe+1,1) = MSC_value;  % Result (not to be ported)
    R.MSC_wind(iframe+1,1)  = MSC_wind ;  % Result (not to be ported)
    
    
    %----------------------------------------------------------------------
    % `` WNS method (signal sub-band centroid)
    %
    % Parameter  f_ssc, th_ssc, alpha (same as MSC)
    % State      XC11, XC22 (time weighted psd's from MSC)
    % Input      X_vpu, X_mic (already clculated in MSC)
    % Results    SSC, wind_decision
    % Output     flag
    %
    %----------------------------------------------------------------------
    SSC_n_1 = 0 ;
    SSC_d_1 = 0 ;
    SSC_n_2 = 0 ;
    SSC_d_2 = 0 ;
    SSC_f = S.SSC_f ;
    % i_SSC_f is index of element at SSC_f
    for mm = 1:i_SSC_f
        SSC_n_1 = SSC_n_1 + XC_11(mm)*mm*fs/nfft;
        SSC_d_1 = SSC_d_1 + XC_11(mm);
        SSC_n_2 = SSC_n_2 + XC_22(mm)*mm*fs/nfft;
        SSC_d_2 = SSC_d_2 + XC_22(mm);
    end
    SSC_1       = SSC_n_1 / SSC_d_1 ;           % [Hz] (scalar)
    SSC_1_value = (SSC_f - SSC_1) / SSC_f ;     % (scalar)
    SSC_2       = SSC_n_2 / SSC_d_2 ;
    SSC_2_value = (SSC_f - SSC_2) / SSC_f ;
    
    % use VPU only
    if (SSC_1_value > S.th_ssc)
        SSC_wind = 1 ;  % classify as wind
    else
        SSC_wind = 0 ;  % not wind
    end
    
    R.SSC_1_value(iframe+1,1) = SSC_1_value ; % Result (not to be ported)
    R.SSC_2_value(iframe+1,1) = SSC_2_value ; % Result (not to be ported)
    R.SSC_wind(iframe+1,1)    = SSC_wind    ; % Result (not to be ported)
    
    
    %----------------------------------------------------------------------
    % RMS diff method (Better: RMS level diff method)
    %----------------------------------------------------------------------
    %rms_vpu    = rms(xframe_vpu);               % (vpu) scalar
    rms_vpu    = sqrt(sum(xframe_vpu.^2)/nwin);               % (vpu) scalar
    rms_dB_vpu = 20*log10(rms_vpu+1.0e-12);     % +1.0e-12 to avoid log(0)!!!
    %rms_mic    = rms(xframe_mic);               % (mic) scalar
    rms_mic    = sqrt(sum(xframe_mic.^2)/nwin);               % (vpu) scalar
    rms_dB_mic = 20*log10(rms_mic+1.0e-12);
    
    rms_diff = S.RMSD_alpha*rms_diff + (1-S.RMSD_alpha)*(rms_dB_mic - rms_dB_vpu) ;
    
    % Results for plotting only
    RMSD_value = rms_diff ;
    % use RMS diff only
    if (RMSD_value > 14.2)      %18.5 for Sigma data, 7.2 for Sage data, 14.5 for Long file
        RMSD_wind = 1 ;  % classify as wind
    else
        RMSD_wind = 0 ;  % not wind
    end
    
    R.RMSD_value(iframe+1,1) = RMSD_value ; % Result (not to be ported)
    R.RMSD_wind(iframe+1,1)  = RMSD_wind  ; % Result (not to be ported)
    
    
    %----------------------------------------------------------------------
    % WNS method selection and application
    %----------------------------------------------------------------------
    
    %     %%%AA = ((1 - MSC_ave) > th_msc);    % for Long file
    %     AA = ((1 - MSC_ave) > 0.88);      % 9/3/2021, 0.92 for OTOG2 case01, 0.90 for case2
    %     BB = (SSC_1_flag > 0.4);
    %     %%%CC = (mean(rms_diff) > 14.5);     % 18.5 for Sigma data, 7.2 for Sage data, 14.5 for Long file
    %
    %     CC = (mean(rms_diff) > 16.5);
    %     %if ~(CC)     % for Long file
    %     %if ~(AA)     % for first 3 cases
    %     if ~(CC)
    %         detect_result((m)*nhop+(1:nhop),1) = 0.5;
    %         AllT_rate(m+1,1) = 0;
    %
    %         frameGain = 1.0;
    %         counter = 0;
    %         %elseif AA || BB || (CC)   % for Long file
    %     else
    %         detect_result((m)*nhop+(1:nhop),1) = 1.0;
    %         AllT_rate(m+1,1) = 1;
    %
    %         counter = counter + 1;
    %         if (counter >= 10)
    %             frameGain = 0.15;    %0.1
    %         end
    %     end
    
    switch S.WNS_method % 'MSC'|'SSC'|'RMSD'|'none'
        case 'MSC'
            wind = MSC_wind ;
        case 'SSC'
            wind = SSC_wind ;
        case 'RMSD'
            wind = RMSD_wind ;
        case 'none'
            wind = 0 ; % (no wind)
        otherwise
            error('Unknown WNSmethod')
    end
    R.wind(iframe+1,1) = wind ;  % Result (not to be ported)
    
    if wind
        WNS_counter = WNS_counter + 1 ;
        if (WNS_counter >= S.WNS_counter_th)
            frameGain = S.WNS_frameAttenuation ;
        else
            frameGain = 1.0 ;
        end
    else
        WNS_counter = 0 ;
        frameGain = 1.0 ;
    end
    
    %----------------------------------------------------------------------
    % WNS output signal processing
    % (Why apply the same frame gain in both channels????)
    %----------------------------------------------------------------------
    %     %y_left = (xframe_l / 20) * frameGain;            % for Long file /20
    %     y_left = (xframe_l / 2.0) * frameGain;
    %     y(m*nhop+(1:M), 2) = y(m*nhop+(1:M), 2) + y_left;   % Ch2 to match with input
    %     %y_right = (xframe_r / 20) * frameGain;           % for Long file /20
    %     y_right = (xframe_r / 2.0) * frameGain;
    %     y(m*nhop+(1:M), 1) = y(m*nhop+(1:M), 1) + y_right;
    
    xframe_vpu = xframe_vpu * frameGain ;
    xframe_mic = xframe_mic * frameGain ;
    
    
    %--------------------------------------------------------------------------
    % Add HPF for acoustic mic (wind noise rejecting first stage in signal chain)
    % *** Make this frame based also!
    % *** FB processing order:  WNS -> HP -> SE
    %--------------------------------------------------------------------------
    if micHPfilterInFrame
        [xframe_mic,HPmic_zf] = filter(HPmic_b, HPmic_a, xframe_mic);
    end


    
    
    
    %     %----------------------------------------------------------------------
    %     % Voice activity detection? (currently not used by the algorithm)
    %     %----------------------------------------------------------------------
    %     % Coherence (with "smoothing", i.e. exponential time weighting)
    %     XC_11 = alpha*XC_11 + (1-alpha)*(X_vpu.*conj(X_vpu));
    %     XC_22 = alpha*XC_22 + (1-alpha)*(X_mic.*conj(X_mic));
    %     XC_12 = alpha*XC_12 + (1-alpha)*(X_vpu.*conj(X_mic));
    %     % MSC (magnitude squared coherence)
    %     MSC = (XC_12.*conj(XC_12)) ./ (XC_11.*XC_22);
    %
    %     % average of lower frequency 0 - 500Hz, 500/(Fs/nfft) -> 21.xx -> 22
    %     MSC_ave = sum(MSC(1:44))/44;
    %     vad_decision((iframe)*nhop+(1:nhop),1) = 1 - MSC_ave;
    %     % set flags
    %     if ((1 - MSC_ave) > th_msc)
    %         detect_result((iframe)*nhop+(1:nhop),1) = 0.5;
    %     else
    %         detect_result((iframe)*nhop+(1:nhop),1) = 1.0;
    %     end
    
    switch S.SE_method % 'STEQ'|'LP'|'HPLP'|'none'
        case 'HPLP'
            %--------------------------------------------------------------
            % Speech enhancement by equalizer
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % STFT with exponential time weighting
            
            % Apply window and fft
            xw_vpu = win .* xframe_vpu;                         % Apply window
            xw_mic = win .* xframe_mic;
            %xwzp_vpu = [xw_vpu(Mo2+1:nwin);zp;xw_vpu(1:Mo2)];      % in zero-phase form!!
            %xwzp_mic = [xw_mic(Mo2+1:nwin);zp;xw_mic(1:Mo2)];
            xwzp_vpu = [xw_vpu;zp];     % to use IFFT directly!!!
            xwzp_mic = [xw_mic;zp];
            X_vpu = fft(xwzp_vpu);
            X_mic = fft(xwzp_mic);
            
            % Exponentially time weighted STFT
            % (Note that setting alpha=0 gives the non-weighted case) 
            Xtw_vpu = S.STEQ_alpha_X*Xtw_vpu + (1-S.STEQ_alpha_X)*X_vpu ;
            Xtw_mic = S.STEQ_alpha_X*Xtw_mic + (1-S.STEQ_alpha_X)*X_mic ;
            
            % Spectral smoothing using time weighted STFT

            % average spectrum of before and after j frequency bins, the first and
            % last jj values are not averaged
            jj = 5;    %3;
            for mm = jj+1:length(Xtw_vpu)-jj
                avgAmpX_vpu(mm) = sum(abs(Xtw_vpu(mm-jj:mm+jj)))/(2*jj+1);
                avgAmpX_mic(mm) = sum(abs(Xtw_mic(mm-jj:mm+jj)))/(2*jj+1);
            end
            avgAmpX_vpu = [abs(Xtw_vpu(1:jj)); avgAmpX_vpu(jj+1:length(Xtw_vpu)-jj); abs(Xtw_vpu(length(Xtw_vpu)-jj+1:length(Xtw_vpu)))];
            avgAmpX_mic = [abs(Xtw_mic(1:jj)); avgAmpX_mic(jj+1:length(Xtw_mic)-jj); abs(Xtw_mic(length(Xtw_mic)-jj+1:length(Xtw_mic)))];
            
            % Amplitude ratio (between frequency smoothed sensor spectra)
            Heq = avgAmpX_mic ./ avgAmpX_vpu;
            
            % Equalized contact mic output
            %equalizedX_vpu = (abs(X_vpu) .* Heq) .* exp(1j .* angle(X_vpu));
            equalizedX_vpu = Heq .* X_vpu ;
            
            % Log Spectral Deviation
            R.LSD(iframe+1, 1) = sqrt(sum((log10(equalizedX_vpu.*conj(equalizedX_vpu))-log10(X_mic.*conj(X_mic))) .^2));
            
            % Processed/enhanced output
            y_enhanced = real(ifft(equalizedX_vpu)) ;
            y_ref = real(ifft(X_mic)) ;
            if 1 % Lin's code
                Ifft = iframe*nhop+(1:nfft) ;
                y(Ifft, 1) = y(Ifft, 1) + y_enhanced;
                y(Ifft, 2) = y(Ifft, 2) + y_ref;
            else % as in WNS code
                y(Iframe, 1) = y(Iframe, 1) + y_enhanced(1:nwin) ;
                y(Iframe, 2) = y(Iframe, 2) + y_ref(1:nwin);
            end
            
        case 'LP'
            %--------------------------------------------------------------
            % Speech enhancement by linear predictor
            %--------------------------------------------------------------
            
            % Apply time window
            xw_vpu = win .* xframe_vpu;
            xw_mic = win .* xframe_mic;
            
            % calculate LPC coefficients
            A_vpu = lpc(xw_vpu, S.LP_order) ;
            A_mic = lpc(xw_mic, S.LP_order) ;

            x_res = filter( A_vpu, A_mic, xw_vpu) ;

            % Processed/enhanced output
            y(Iframe, 1) = y(Iframe, 1) + x_res ;
            y(Iframe, 2) = y(Iframe, 2) + xw_mic ; % from Lin's code (use time windowed frames????)
            
        case 'HPLP'
            %--------------------------------------------------------------
            % Speech enhancement by cross-over filter
            %--------------------------------------------------------------
            
            % Apply time window
            xw_vpu = win .* xframe_vpu ;
            xw_mic = win .* xframe_mic ;
            
            % HP for acoustic mic; LP for contact mic
            x_lp = filter(Bl, Al, xw_vpu);    % CH: filtering without saving internal status of filter for next frame!!
            x_hp = filter(Bh, Ah, xw_mic);
            % Combine LP and HP signals
            x_com = x_lp + x_hp;

            % Processed/enhanced output
            y(Iframe, 1) = y(Iframe, 1) + x_com ;
            y(Iframe, 2) = y(Iframe, 2) + xw_mic ; % from Lin's code (use time windowed frames????)
            
        case 'none'
            y(Iframe, 1) = xframe_vpu ; % is this correct???
            y(Iframe, 2) = xframe_mic ;
            
        otherwise
            error('Non-valid SE method.')
    end
        
    if rem(iframe,round(fs/nhop))==0
        elapsedTime = toc ;
        estimatedRemainingTime = elapsedTime/iframe*(nframes-iframe) ;
        waitbar(iframe/nframes,fProgress,['Estimated remaining time: ' num2str(round(estimatedRemainingTime)) ' s.']) ;
    end
end
elapsedTime = toc ;
close(fProgress)
disp(['Processing time/real time: ' num2str(elapsedTime/(nx/fs))])

y = y(1:nx,:) ; % Make input and output signals equal length
R.fs            = fs ;

Sout = S ;    % NOTE Any non-valid field values of S may have been changed

end

