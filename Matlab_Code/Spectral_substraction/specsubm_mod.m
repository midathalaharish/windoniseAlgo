function [ss,po]=specsubm_mod(s,fs,s2,p)
%SPECSUBM performs speech enhancement using spectral subtraction [SS,PO]=(S,FS,P,S2)
%
% S is the input mixture, speech + noise
% S2 is the clean speech (modification of original SPECSUBM)
% FS is the sampling rate of S and S2
% P is the parameter vector with details noted below:
%
% implementation of spectral subtraction algorithm by R Martin (rather slow)
% algorithm parameters: t in seconds, f in Hz, k* dimensionless
% 1: tg = smoothing time constant for signal power estimate (0.04): high=reverberant, low=musical
% 2: ta = smoothing time constant for signal power estimate
%        used in noise estimation (0.1)
% 3: tw = fft window length (will be rounded up to 2^nw samples)
% 4: tm = length of minimum filter (1.5): high=slow response to noise increase, low=distortion
% 5: to = time constant for oversubtraction factor (0.08)
% 6: fo = oversubtraction corner frequency (800): high=distortion, low=musical
% 7: km = number of minimisation buffers to use (4): high=waste memory, low=noise modulation
% 8: ks = oversampling constant (4)
% 9: kn = noise estimate compensation (1.5)
% 10:kf = subtraction floor (0.02): high=noisy, low=musical 0.01 < sub
% floor < 0.05
% 11:ko = oversubtraction scale factor (4): high=distortion, low=musical
%
% Refs:
%    (a) R. Martin. Spectral subtraction based on minimum statistics. In Proc EUSIPCO, pages 1182-1185, Edinburgh, Sept 1994.
%    (b) R. Martin. Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
%
%      Copyright (C) Mike Brookes 2004
%      Version: $Id: specsubm.m,v 1.4 2007/05/04 07:01:39 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if nargin<3
  po=[0.04 0.1 0.032 1.5 0.08 400 4 1 1.5 0.02 4].'; 
  s2=s;
else po=p; 
end
 
 
ns=length(s);
s2 = s2(1:ns);
ts=1/fs;
ss=zeros(ns,1);
 
ni=1024;%pow2(nextpow2(fs*po(3)/po(8))); % po(3): fft window length
ti=ni/fs;
nw=ni*po(8); % po(8): oversampling constant
nf=1+floor((ns-nw)/ni);
nm=ceil(fs*po(4)/(ni*po(7))); % po(4): length min filter, po(7) minBufferCount
 
%win=0.5*hamming(nw+1)/1.08;win(end)=[];
 
%win=sqrt(hamming(nw+1)); win(end)=[]; % for now always use sqrt hamming window
win=(hann(nw)); % for now always use sqrt hamming window
##win=win/sqrt(sum(win(1:ni).^2));       % normalize to give overall gain of 1
 
zg=exp(-ti/po(1)); % po(1) smoothing time constant
za=exp(-ti/po(2)); % po(2) smoothing time constant in noise estimate
zo=exp(-ti/po(5)); % po(5) time constant for oversubtraction factor 
 
px=zeros(nw/2,1);
pxs = px;
pxn=px;
os=px;
pn=px;
mb=ones(nw/2,po(7))*nw/2; % po(7) minBufferCount
im=0;
osf=po(11)*((1:nw/2).'*fs/(nw*po(6))).^(-1); % po(6) oversubtraction corner frequency, po(11) oversubtraction scale factor
 
imidx=[13 21]';
x2im=zeros(length(imidx),nf);
osim=x2im;
pnim=x2im;
pxnim=x2im;
qim=x2im;
 
q = px;
 
for is=1:nf
   idx=(1:nw)+(is-1)*ni;
 
   x=fft(s(idx).*win);
   xs=fft(s2(idx).*win); 
   x2=x.*conj(x);
   xs2=xs.*conj(xs);
   
   pxn=za*pxn+(1-za)*x2(1:nw/2,1); % signal power estimate
   im=rem(im+1,nm);
   if im % always if im ~= 0
      mb(:,1)=min(mb(:,1),pxn);
   else % only if im == 0 (add to MinBuffers)
      mb=[pxn,mb(:,1:po(7)-1)]; % po(7) minBufferCount
   end
   pn=po(9)*min(mb,[],2); % po(9) noise estimate compensation
   %os : oversubtraction factor
   os=zo*os+(1-zo)*(1+osf.*pn./(pn+pxn));
    
   px= zg*px+(1-zg)*x2(1:nw/2,1);
   pxs= zg*pxs+(1-zg)*xs2(1:nw/2,1);
   
   if any(x2 == 0) | any(pn < 0) | any(px == 0) | any(os.*pn < 0)
       q = 0;
   else
       q=max(po(10)*sqrt(pn), 1-sqrt(os.*pn./px)); % po(10) subtraction floor
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
##   D = ifft(xs.*abs(nw_q));
   D = ifft(xs).*(sum(nw_q)/length(nw_q)));

 
%   Y= xs.*(1-sqrt(os.*pn./xs2));
%   Y(Y<sqrt(po(10)*pn))=sqrt(po(10)*pn);
%    D = ifft(Y)
;
  
   ss(idx)=ss(idx)+D.*win;  
##     plot(s(idx),'b');hold on;plot(real(ss(idx)),'r');
##     hold off;
end
 
##ss = ss./max(ss)*0.9;
 
if nargout==0
   soundsc([s; ss],fs);
end