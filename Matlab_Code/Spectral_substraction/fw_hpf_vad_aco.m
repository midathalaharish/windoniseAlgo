function [Y,R,Paco] = fw_hpf_vad_aco(X,fs,nxFrame,P)
%[Y,R] = fw_hpf_vad_hplp(X,fs)

% === Set up frame handling ===============
nX = size(X,1) ; % [Sa] length of input signal(s)
Y = zeros(nX,1) ; % (assume single output)
nxFrame =2048
% === Initialize processing blocks ========


Shpf            = [] ; % State struct must be initialized as empty
Shsf            = [] ; % State struct must be initialized as empty
Svad            = [] ; % State struct must be initialized as empty
Sss            = [] ; % State struct must be initialized as empty
R.vad.VAD       = zeros(nX,1) ;
R.vad.E         = zeros(nX,1) ;
Saco            = [] ; % State struct must be initialized as empty
R.aco.NE_dB     = zeros(nX,1) ;
R.aco.Avg_NE_dB = zeros(nX,1) ;
R.aco.ixf       = zeros(nX,1) ;
R.ss            = zeros(nX,1) ;

if nargin<4

% --- hpf contact mic --------------------
  P.tLookahead = 0 ;
  P.nOverlap   = 0 ; % if overlap is used there is a need for an output history state variable
  P.f0         = 80 ; % [Hz] cut-off frequency
  P.n          = 6 ;  % filter order
  
  
  % --- vad --------------------
  
  P.tLookahead    = 0 ;
  P.nOverlap      = 0 ; % if overlap is used there is a need for an output history state variable
  P.bpf_n         = 2 ;
  P.bpf_f1        = 100 ;   % [Hz]
  P.bpf_f2        = 1500 ;  % [Hz]
  P.tA            = 0.001 ; % time constant [s] for time weighted power
  P.vad_threshold = 1e-4 ;
  P.tHang         = 0.1 ;   % hangover time [s]
  
  % --- high shelving filter --------------------
  P.tLookahead = 0 ;
  P.nOverlap   = 0 ; % if overlap is used there is a need for an output history state variable
##  [b_hs, a_hs] = biquad_filter('highShelf', -12, 700, fs, 0.7);  % high shelving filter for contact mic
##  P.b_hs         = b_hs;
##  P.a_hs         = a_hs;
##  P.n_hs         = 2;
  
  % --- adaptive crossover --------------------
  
  P.tLookahead   = 0 ;
  P.nOverlap     = 0 ; % if overlap is used there is a need for an output history state variable
  P.tB           = 0.001 ;                             % [s] time constant for ambient noise power estimation
  P.Xc           = [ 80, 315, 630, 1250, 2500, 5000] ; % [Hz] frequencies to define crossover
  P.Xc_thres_dB  = [-60 -55 -50 -45 -40 -35] ;         % [dBFS] crossover frequency selection thresholds
  P.n            = 4 ;                                 % x-over filter order

  % --- spectral subtraction --------------------
	P.tLookahead   = 0.010667 ;
  P.nOverlap     = 512 ; % if overlap is used there is a need for an output history state variable
  P.tg = 0.04 ;	%smoothing time constant for signal power estimate (0.04): high=reverberant, low=musical
	P.ta = 0.1 ;   	%smoothing time constant for signal power estimate, used in noise estimation (0.1)
	P.tw = 0.032 ;  %fft window length (will be rounded up to 2^nw samples)
	P.tm = 1.5 ;    %length of minimum filter (1.5): high=slow response to noise increase, low=distortion
	P.to = 0.08 ;	%time constant for oversubtraction factor (0.08)
	P.fo = 400 ;	%oversubtraction corner frequency (800): high=distortion, low=musical
	P.km = 4;		%number of minimisation buffers to use (4): high=waste memory, low=noise modulation
	P.ks = 4 ;		%oversampling constant (4)
	P.kn = 1.5 ;	%noise estimate compensation (1.5)
	P.kf = 0.02 ;	%subtraction floor (0.02): high=noisy, low=musical 0.01 < sub; % floor < 0.05
	P.ko = 4 ;		%oversubtraction scale factor (4): high=distortion, low=musical
	P.ni = pow2(nextpow2(fs*P.tw/P.ks)); % po(3): fft window length
	P.nw = P.ni*P.ks; % po(8): oversampling constant
end

Phpf = P ;
Pvad = P ;
Phsf = P ;
Paco = P ;
Pss  = P ;
% =========================================

% === Frame processing ====================

ip = 1 ; % pointer to first sample in frame
dataAvailable = (ip + nxFrame -1) < nX ;

% --- Set up progress bar ---
nLoop = round(nX/nxFrame) ;
iLoop = 0 ;
fProgress = waitbar(0,['Estimated remaining time: -' ' s.']) ;


% --- Loop -----------------------------------
tic
while dataAvailable
    
    % --- Get new data -----------------------
    Inew = ip+(1:nxFrame)-1 ;
    x_contact = X(Inew,1) ;
    x_acoustic = X(Inew,1) ;
    
		if(0)
    % --- hpf contact mic -------
    [x_contact,Shpf,~] = fb_hpf(x_contact,fs,Shpf,Phpf) ;
    % Store states for graphing and analysis
    %R.fb_hpf.z(Inew,1)  = Rz.z(1:nFrame) ;
    
    % --- vad --------------------
    [VAD,Svad,Rvad] = fb_vad(x_contact,fs,Svad,Pvad) ;
    % Store states for graphing and analysis
    R.vad.VAD(Inew,1)  = VAD(1:nxFrame) ;
    R.vad.E(Inew,1)    = Rvad.E(1:nxFrame) ;
    
    %include EQ here (high shelving for contact mic)
    [x_contact_hs,Shsf,~] = fb_hsf(x_contact,fs,Shsf,Phsf);
    
   
    % --- adaptive crossover (hplp) --------------------
    [x_mix,Saco,Raco, Paco] = fb_adaptive_crossover([x_contact_hs,VAD,x_acoustic],fs,Saco,Paco) ;
    % Store states for graphing and analysis
    R.aco.NE_dB(Inew,1)     = Raco.NE_dB ;
    R.aco.Avg_NE_dB(Inew,1) = Raco.Avg_NE_dB ;
    R.aco.ixf(Inew,1)       = Raco.ixf ;
		end
		% --- spectral subtraction --------------------
%		[y,Sns,Rns,Pns] = frb_ns([x_contact_hs,VAD,x_acoustic],fs,Sns,Pns) ;
		
		[y,Sss,Rss]=frb_ns([x_contact,x_acoustic],fs,Sss) ;
%		[y,S,R]=frb_ns(xFrame,fs,S,P)	
%		[y,Sss,Rss] = frb_SS([x_contact,x_acoustic],fs,Sss,Pss) ;
		
%    y = x_mix ;
    
    
    % --- Output audio data ------------------
    Y(Inew) = y(1:nxFrame) ;
    
    % ----------------------------------------
    % Check if more data is available from main input
    dataAvailable = (ip + 2*nxFrame -1) < nX ;
    
    % Move pointer to beginning of next frame
    ip = ip + nxFrame ;
    
    % --- Update progress bar ---
    iLoop = iLoop+1 ;
    if rem(iLoop,round(fs/(2*nxFrame)))==0
        elapsedTime = toc ;
        estimatedRemainingTime = elapsedTime/iLoop*(nLoop-iLoop) ;
        waitbar(iLoop/nLoop,fProgress,['Estimated remaining time: ' num2str(round(estimatedRemainingTime)) ' s.']) ;
    end
    % ----------------------------
    
end
% --- Clean up progress bar ---
elapsedTime = toc ;
close(fProgress)
disp(['Processing time/real time: ' num2str(elapsedTime/(nX/fs))])
% -----------------------------

end