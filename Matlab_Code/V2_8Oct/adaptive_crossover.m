function [x_com, xw_l, xw_r, noise_est_frame,Wn] = adaptive_crossover(xframe_l, xframe_r, Xc, Fs, noise_est_var, vad_signal, M, Win)

alpha_exp = 1/(Fs*0.01);

%% HP & LP coeffcients

[Bl1, Al1]  = butter(2, Xc(1)/(Fs/2), 'low')
[Bh1, Ah1]  = butter(2, Xc(1)/(Fs/2), 'high')

[Bl2, Al2]  = butter(2, Xc(2)/(Fs/2), 'low')
[Bh2, Ah2]  = butter(2, Xc(2)/(Fs/2), 'high')

[Bl3, Al3]  = butter(2, Xc(3)/(Fs/2), 'low')
[Bh3, Ah3]  = butter(2, Xc(3)/(Fs/2), 'high')

[Bl4, Al4]  = butter(2, Xc(4)/(Fs/2), 'low')
[Bh4, Ah4]  = butter(2, Xc(4)/(Fs/2), 'high')

[Bl5, Al5]  = butter(2, Xc(5)/(Fs/2), 'low')
[Bh5, Ah5]  = butter(2, Xc(5)/(Fs/2), 'high')

[Bl6, Al6]  = butter(2, Xc(6)/(Fs/2), 'low')
[Bh6, Ah6]  = butter(2, Xc(6)/(Fs/2), 'high')

%% estimate noise when VAD is 0

q = noise_est_var;

for k = 1:M
  if(vad_signal == 0)
  q  =  (1-alpha_exp)*q + alpha_exp*xframe_r(k).^2;
endif
noise_est_frame(k)  =   q;
endfor

noise_est_var_dB  =   10*log10(q + 1.0e-12); % to avoid log 0

if (noise_est_var_dB  <= -60)
  Wn = Xc(1);
  Bl = Bl1; Al = Al1; Bh = Bh1; Ah = Ah1;
  elseif (noise_est_var_dB  <= -50)
  Wn = Xc(2);
  Bl = Bl2; Al = Al2; Bh = Bh2; Ah = Ah2;
  elseif (noise_est_var_dB  <= -40)
  Wn = Xc(3);
  Bl = Bl3; Al = Al3; Bh = Bh3; Ah = Ah3;
  elseif (noise_est_var_dB  <= -30)
  Wn = Xc(4);
  Bl = Bl4; Al = Al4; Bh = Bh4; Ah = Ah4;
  elseif (noise_est_var_dB  <= -20)
  Wn = Xc(5);
  Bl = Bl5; Al = Al5; Bh = Bh5; Ah = Ah5;
else 
  Wn = Xc(6);
  Bl = Bl6; Al = Al6; Bh = Bh6; Ah = Ah6;
endif

%% windowing

xw_l  =   Win.* xframe_l;
xw_r  =   Win.* xframe_r;  

%% filtering

x_lp  =   filter(Bl, Al, xw_l);
x_hp  =   filter(Bh, Ah, xw_r);

% add both

x_com   =  x_lp +  x_hp;
 
end
