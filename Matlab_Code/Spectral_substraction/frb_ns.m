%function [ss,po]=ss_ns_mod1(s,fs,s2,p)
function [y,S,R]=frb_ns(xFrame,fs,S,P)	
 
%[y,S,R] = frb_hpf(xFrame,fs,S,P)
%
%    xFrame   new audio samples 
%    fs       sample rate
%    S        state struct
%    P        parameter struct
%    y        block output (audio samples)
%    R        results struct (for visualizing internal quantities) 

% implementation of spectral subtraction algorithm by R Martin




nxFrame = size(xFrame,1) ; % number of new audio samples
nCh     = size(xFrame,2) ; % number of input channels

if nargin<4
    P = [] ;
end
% Default parameter value set
if isempty(P)
	P.tLookahead = 0 ;
  P.nOverlap = 1536-nxFrame  ; % if overlap is used there is a need for an output history state variable
	
	P.tg = 0.04 ;	%smoothing time constant for signal power estimate (0.04): high=reverberant, low=musical
	P.ta = 0.1 ;   	%smoothing time constant for signal power estimate, used in noise estimation (0.1)
	P.tw = 0.032 ;  %fft window length (will be rounded up to 2^nw samples)
	P.tm = 1.5 ;    %length of minimum filter (1.5): high=slow response to noise increase, low=distortion
	P.to = 0.08 ;	%time constant for oversubtraction factor (0.08)
	P.fo = 400 ;	%oversubtraction corner frequency (800): high=distortion, low=musical
	P.km = 4;		%number of minimisation buffers to use (4): high=waste memory, low=noise modulation
	P.ks = 1 ;		%oversampling constant (4)
	P.kn = 1.5 ;	%noise estimate compensation (1.5)
	P.kf = 0.02 ;	%subtraction floor (0.02): high=noisy, low=musical 0.01 < sub; % floor < 0.05
	P.ko = 4 ;		%oversubtraction scale factor (4): high=distortion, low=musical
	
	P.ni = 1536;%pow2(nextpow2(fs*P.tw/P.ks)); % po(3): fft window length
	P.nw = P.ni*P.ks; % po(8): oversampling constant
	
end
nLookahead 	= round(fs*P.tLookahead) ;
nHistory 	= P.nOverlap + nLookahead  ;

% Initial block state
if isempty(S)
		S.pxn  = zeros(P.nw/2,1) ;
		S.os  = zeros(P.nw/2,1) ;
		S.px  = zeros(P.nw/2,1) ;
		S.pxs  = zeros(P.nw/2,1) ;
		
end

% Add new audio samples and update input history
x =xFrame(:,1) ;

for is=1:2					
%--------------------------------------------------------------------------
% Process audio data
idx = (1:1536) + (is -1) *512;
s  = x(idx) ;
s2 = s ;
ns=length(s);
s2 = s2(1:ns);
ts=1/fs;
ss=zeros(ns,1);
 
%ni=pow2(nextpow2(fs*P.tw/P.ks)); % po(3): fft window length
%ni = ni/2;
ti=P.ni/fs;
%nw=ni*P.ks; % po(8): oversampling constant
nf=1+floor((ns-P.nw)/P.ni);
nm=ceil(fs*P.tm/(P.ni*P.km)); % po(4): length min filter, po(7) minBufferCount
 
%win=0.5*hamming(nw+1)/1.08;win(end)=[];
 
%win=sqrt(hamming(nw+1)); win(end)=[]; % for now always use sqrt hamming window
win=sqrt(hann(P.nw)); % for now always use sqrt hamming window
win=win/sqrt(sum(win(1:P.ni).^2));       % normalize to give overall gain of 1
 
zg=exp(-ti/P.tg); % po(1) smoothing time constant
za=exp(-ti/P.ta); % po(2) smoothing time constant in noise estimate
zo=exp(-ti/P.to); % po(5) time constant for oversubtraction factor 
 
%px=zeros(P.nw/2,1);
%pxs = px;
%pxn=px;
%os=px;
pn=S.px;
mb=ones(P.nw/2,P.km)*P.nw/2; % po(7) minBufferCount
im=0;
osf=P.ko*((1:P.nw/2).'*fs/(P.nw*P.fo)).^(-1); % po(6) oversubtraction corner frequency, po(11) oversubtraction scale factor
%osf=po(11)*((1:nw/2).'*fs/(nw*po(6))).^(-1); 
%imidx=[13 21]';
%x2im=zeros(length(imidx),nf);
%osim=x2im;
%pnim=x2im;
%pxnim=x2im;
%qim=x2im;
 
q = S.px;
 
%for is=1:nf
%   idx=(1:nw)+(is-1)*ni;
 
   x=fft(s.*win);
   xs=fft(s2.*win); 
%	 x = x(1:(floor(length(x)/2)+1));
%	 xs = xs(1:(floor(length(xs)/2)+1));
   x2=x.*conj(x);
   xs2=xs.*conj(xs);
   
   S.pxn=za*S.pxn+(1-za)*x2(1:P.nw/2,1); % signal power estimate
   im=rem(im+1,nm);
   if im % always if im ~= 0 fix this
      mb(:,1)=min(mb(:,1),S.pxn);
   else % only if im == 0 (add to MinBuffers)
      mb=[S.pxn,mb(:,1:P.km-1)]; % po(7) minBufferCount
   end
   pn=P.kn*min(mb,[],2); % po(9) noise estimate compensation
   %os : oversubtraction factor
   S.os=zo*S.os+(1-zo)*(1+osf.*pn./(pn+S.pxn));
    
   S.px= zg*S.px+(1-zg)*x2(1:P.nw/2,1);
   S.pxs= zg*S.pxs+(1-zg)*xs2(1:P.nw/2,1);
   
   if any(x2 == 0) | any(pn < 0) | any(S.px == 0) | any(S.os.*pn < 0)
       q = 0;
   else
       q=max(P.kf*sqrt(pn), 1-sqrt(S.os.*pn./S.px)); % po(10) subtraction floor
%        q(1)=0; 
%         q=filter(1/3*[1 1 1],1,q);%q=filter(0.8,[1 -0.2],q);
%        q=filter(1/3*[1 1 1],1,q(end:-1:1));
%        q=q(end:-1:1);
%        q(1)=0;
%        plot(q);drawnow;
   end   
   %D = ifft(x.*q);
   
   nw_q = flip(q);
   nw_q = repelem(nw_q,2);
   D = ifft(xs.*abs(nw_q));
% xss = xs.*abs(q);
%	 even = true;
%	 n11 = 0; % the output length
%   s11 = 0; % the variable that will hold the index of the highest
%	 % frequency below N/2, s = floor((n+1)/2)
%	 if (even)
%	 n11 = 2*(length(xss)-1);
%	 s11 = length(xss)-1;§
%	 else
%	 n11 = 2 * (length(xss) - 1 )+1;
%	 s11 = length(xss);
%	 endif
%	 xn = zeros(1,n11);
%	 xn(1:length(xss)) = xss;
%	 xn(length(xss):n11) = conj(xss(s11:-1:1));
%	D = ifft(xn);
  
##  ss=ss+real(D).*win;  
	##ss = ss./max(ss)*0.9;
%end
	 y = ss;
##     plot(s(idx),'b');hold on;plot(real(ss(idx)),'r');
##     hold off;

%--------------------------------------------------------------------------
ss(idx)=real(ss(idx)); 

% Result output
%R.nLookahead = nLookahead ;

% Audio output
end
y=ss;
end
 