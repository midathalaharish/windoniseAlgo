function [y,S,R]=frb_SS(xFrame,fs,S,P)										 
%[y,S,R] = frb_hpf(xFrame,fs,S,P)
%
%    xFrame   new audio samples 
%    fs       sample rate
%    S        state struct
%    P        parameter struct
%    y        block output (audio samples)
%    R        results struct (for visualizing internal quantities) 

% implementation of spectral subtraction algorithm by R Martin

% 2021Nov17 spectral subtraction based on minimum statistics
% 2021Nov16 Naming overhaul: "fb_" changed to "frb_" and "fw_" to "frw_"
% 2021Nov03rd first freeze
%   Copyright 2021 Sigma Connectivity AB ^

if nargin<4
    P = [] ;
end
% Default parameter value set
if isempty(P)
	P.tLookahead = 0 ;
  P.nOverlap = 0 ; % if overlap is used there is a need for an output history state variable
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
end
nLookahead 	= round(fs*P.tLookahead) ;
nHistory 	= P.nOverlap + nLookahead  ;

% Initial block state
if isempty(S)
    S.xHistory = zeros(nHistory,1) ;
    %S.z = zeros(P.n,1) ; %  initial/final conditions z of the filter delays
end

% Add new audio samples and update history
nxFrame = length(xFrame) ; % number of new audio samples
x = [S.xHistory; xFrame] ;
S.xHistory = x(nxFrame + (1:nHistory)) ;

					
%--------------------------------------------------------------------------
% Process audio data
x_contact  = x(:,1) ;
%vad        = x(:,2) ;
%x_acoustic = x(:,3) ;
%fft_length = fs*(P.tw/(P.ks-1));									% pow2(nextpow2(fs*P.tw/P.ks)); % P.tw: fft window length
fft_length = 256; 															  % changed to match xFrame length
time_vec   = fft_length/fs;
win_len    = fft_length*P.ks/4;										  % oversampling constant
%win_len		 = fft_length;													% changed to match xFrame length
im		     = 0;
nm		     = ceil(fs*P.tm/(fft_length*P.km)); 						% P.tm: length min filter, P.km minBufferCount
window	   = hann(win_len); 									% for now always use sqrt hamming window
%window 	   = window/sqrt(sum(window(1:fft_length:win_len).^2)); 	% normalize to give overall gain of 1
alpha_smo  = exp(-time_vec/P.tg); 										% P.tg smoothing time constant
alpha_ne   = exp(-time_vec/P.ta); 										% P.ta smoothing time constant in noise estimate
alpha_osub = exp(-time_vec/P.to); 										% P.to time constant for oversubtraction factor
%osub	     = P.ko*(1+(0:win_len/2).'*fs/(win_len*P.fo)).^(-1); 		% P.fo oversubtraction corner frequency, P.ko oversubtraction scale factor
osub 			 = P.ko*((1:win_len).'*fs/(win_len*P.fo)).^(-1);		% changed to match length
%min_buf	   = ones(1+win_len/2,P.km)*win_len/2; 						% po(7) minBufferCount
min_buf    = ones(win_len,P.km)*win_len/2;
%mb=ones(1+nw/2,po(7))*nw/2;
ss_out 	   = zeros(nxFrame,1) ;
sp_est 	   = zeros(nxFrame,1) ;
sp_est_smo = zeros(nxFrame,1) ;
ne_est 	   = zeros(nxFrame,1) ;
osub_est   = zeros(nxFrame,1) ;
q          = zeros(nxFrame,1) ;


%for ix = 1:nxFrame
%	x_contact_fft = fft(x_contact.*window);
x_contact_fft = fft(x_contact.*window);
x_contact_fftabs = x_contact_fft.*conj(x_contact_fft);
	
	sp_est = alpha_ne*sp_est+(1-alpha_ne)*x_contact_fftabs; % signal power estimate
	im = rem(im+1,nm);
	if im % always if im ~= 0
		min_buf(:,1) = min(min_buf(:,1),sp_est);
	else % only if im == 0 (add to MinBuffers)
		min_buf = [sp_est, min_buf(:,1:P.km-1)]; % P.km minBufferCount
	end
	ne_est = P.kn*min(min_buf,[],2); % P.kn noise estimate compensation
	osub_est = alpha_osub*osub_est+(1-alpha_osub)*(1+osub.*ne_est./(ne_est+sp_est));
	
	sp_est_smo= alpha_smo*sp_est_smo+(1-alpha_smo)*x_contact_fftabs;
	
	if any(x_contact_fftabs == 0) | any(ne_est < 0) | any(sp_est_smo == 0) | any(osub_est.*ne_est < 0)
		q = 0;
	else
		q = max(P.kf*sqrt(ne_est), 1-sqrt(osub_est.*ne_est./sp_est_smo)); % P.kf subtraction floor
	end
%	q1 = flip(q);
%	q = [q;q1];
%	q = q(1: length(q)-2,1);

%	D11 = x_contact_fft.*abs(q);
%	n11 = 0;
%	s11 = 0;
%	n11 = 2*(length(D11)-1);
%	s11 = length(D11)-1;
%	xn1 = zeros(1, n11);
%	xn1(1:length(D11):n11) = conj(D11(s11:-1:1));
%end;
%q = 1;
	xss = x_contact_fft.*q;
%	D1 = fft_length*ifft(ifftshift(xss),fft_length);		%N*ifft(ifftshift(X),N); 
	D1 = ifft(xss);
	ss_out = ss_out + D1;
%	y = ss_out + D;
%end

%y = ss_out./max(ss_out)*0.1;
y = ss_out;
%y = x_contact;
%--------------------------------------------------------------------------
Iout = nHistory + (1:nxFrame) - nLookahead ;

% Result output
%R.z = z(Iout) ;
R = [] ;

% Audio output
y = y(Iout) ;

end