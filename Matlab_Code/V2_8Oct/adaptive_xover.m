%%%%% ADAPTIVE CROSSOVER FREQUENCIES OF HPF-LPF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Read input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xin, Fs] = audioread (input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    pre treat input with HPF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a]            = butter(6, 80/*Fs/2), 'high');
x_right_acoustic = filter(b,a, xin(:,1))
x_left_vpu       = filter(b,a, xin(:,2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    pre treat vpu signal with BPF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band_freqs = [100,1500];
[b_high, a_high] = butter(2, band_freqs(1)/Fs, 'high');
[b_low, a_low] = butter(2, band_freqs(2)/Fs, 'low');
b2 = conv(b_high, a_low);
a2 = conv(a_high, b_low);
x_vpu_filtered = filter(b2, a2, x_in(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    variables declaration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nooverlap = 256;
window = hann(512, 'periodic');
M = length(window);
w = window(:);
nhop  = M-nooverlap;
nx = length(x_left_vpu);
nframes = 1+fix((nx-M)/nhop);
xframe_l = zeros(M,1);
xframe_r = zeros(M,1);
noise_est_var = 0;
noise_est_frame = zeros(M,1);
noise_est_full = zeros(nframes*nhop+M,1);
filtered_out = zeros(nframes*nhop+M,2);
vad_signal = zeros(M,1);
vad_signal_full = zeros(nframes*nhop+M,1);
vad_threshold = 2e-5;
vad_average = 0;
vad_average_frame   = zeros(M,1); %E_out
vad_average_full = zeros(nframes*nhop+M,1); % E_full
yout = zeros(nframes*nhop+M,2);
t = (0:length(nframes)-1)*1/Fs;
Xc = [80,315,630,1250,2500,5000];



alpha = 1/(Fs*0.01);
vad_ths = 1;

for m = 0:nframes-1
  xframe_l =  x_left_vpu(m*nhop+(1:M));
  xframe_l_filtered = x_vpu_filtered(m*nhop+(1:M));
  xframe_r =  x_right_acoustic(m*nhop+(1:M)); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    VAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [vad_signal, vad_average_frame] = VAD(xframe_l_filtered, vad_average, alpha, M, nhop, vad_threshold );
  vad_average = vad_average_frame(nhop);
  vad_average_full(m*nhop+(1:nhop),1) = vad_average_frame;
  vad_signal_full(m*nhop+(1:nhop),1) = vad_signal;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ADAPTIVE cross over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [x_com, xw_l, xw_r, noise_est_frame] = adaptive_crossover(xframe_l, xframe_r, Xc, Fs, noise_est_var, vad_signal, M, w);
  noise_est_var = noise_est_frame(nhop);
  noise_est_full(m*nhop+(1:M),1)  = noise_est_frame;
  
  yout(m*nhop+(1:M),1) = yout(m*nhop+(1:M),1)+ x_com;
  yout(m*nhop+(1:M),2) = yout(m*nhop+(1:M),2)+ xw_r;
  endfor
  
endfor
