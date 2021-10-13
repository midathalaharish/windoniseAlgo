  %%%%%%%%%%   LPF-HPF method uding adaptive crossover   %%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Description: A method to choose crossover frequencies in HPF-LPF method adaptively.
  %
  %
  % Date:        9/30/2021, initial, based on WNS_Algorithm.m & SE_Algorithm_HP_LP_filter
  %
  % Notes:       9/30/2021, need more fine tune and improvement. Check with
  %                        the paper and other description. More test and
  %                        comparisons
  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %-------------------------------------------------------------------------
  % load input data
  %-------------------------------------------------------------------------
  
  function [x_com, xw_l, xw_r, NE, sum_ne, Wn]  =  adaptive_crossover(xframe_l, xframe_r, Xc, Fs, exp_avg_r, VAD_signal, M, w, sum_ne)

  %--------------------------------------------------------------------------
  % start algorithm
  %--------------------------------------------------------------------------
  tB = 0.01;
  alpha_exponential = 1/(Fs*tB);  % need to be tuned
%  alpha_exponential = 1;  % need to be tuned
%  alpha_exponential = 0;  % need to be tuned
%  sum_ne = 0;
  %--------------------------------------------------------------------------
  % HP and LP Butterworth filters
  %--------------------------------------------------------------------------
  
%  Xc = [ 80, 315, 630, 1250, 2500, 5000];        % frequencies to define crossover
  
  [Bl1,Al1] = butter(4,Xc(1)/(Fs/2),'low');       % Butterworth filter low pass
  [Bh1,Ah1] = butter(4,Xc(1)/(Fs/2),'high');      % Butterworth filter high pass
  
  [Bl2,Al2] = butter(4,Xc(2)/(Fs/2),'low');       % Butterworth filter low pass
  [Bh2,Ah2] = butter(4,Xc(2)/(Fs/2),'high');      % Butterworth filter high pass
  
  [Bl3,Al3] = butter(4,Xc(3)/(Fs/2),'low');       % Butterworth filter low pass
  [Bh3,Ah3] = butter(4,Xc(3)/(Fs/2),'high');      % Butterworth filter high pass
  
  [Bl4,Al4] = butter(4,Xc(4)/(Fs/2),'low');       % Butterworth filter low pass
  [Bh4,Ah4] = butter(4,Xc(4)/(Fs/2),'high');      % Butterworth filter high pass
  
  [Bl5,Al5] = butter(4,Xc(5)/(Fs/2),'low');       % Butterworth filter low pass
  [Bh5,Ah5] = butter(4,Xc(5)/(Fs/2),'high');      % Butterworth filter high pass
  
  [Bl6,Al6] = butter(4,Xc(6)/(Fs/2),'low');       % Butterworth filter low pass
  [Bh6,Ah6] = butter(4,Xc(6)/(Fs/2),'high');      % Butterworth filter high pass
  

   % implement noise estimations in samples
   
%   for k = 1:M
%     if (VAD_signal == 0)                 % estimate the average of acoustic mic in noisy conditions only
%       exp_avg_r   = (1-alpha_exponential)*exp_avg_r + alpha_exponential*(xframe_r(k).^2);
%     else 
%       exp_avg_r   =  exp_avg_r;
%     endif
%     NE(k)         = exp_avg_r;
%    end 
    
    q = exp_avg_r;                        % start with nhop(256) and read the last updated value
    for k = 1:M
      if (VAD_signal == 0)   
    % Compute sample energy
      q = (1-alpha_exponential)*q + alpha_exponential*xframe_r(k).^2;
      endif
    NE(k) = q;
    sum_ne(k) += NE(k);
    end
    
    %% calculate energy  for averaged acoustic mic signal
   
   xframe_r_avg_dB    = 10*log10(q+1.0e-12); % to avoid log(0)   % rename it to nosie eastimatte
   
    
    if(xframe_r_avg_dB <= -60)
     Wn = Xc(1);
      Bl = Bl1; Al = Al1; Bh = Bh1; Ah = Ah1;
      elseif(xframe_r_avg_dB <= -50)
      Wn = Xc(2);
      Bl = Bl2; Al = Al2; Bh = Bh2; Ah = Ah2;
      elseif(xframe_r_avg_dB <= -40)
      Wn = Xc(3);
      Bl = Bl3; Al = Al3; Bh = Bh3; Ah = Ah3;
      elseif(xframe_r_avg_dB <= -30)
      Wn = Xc(4);
      Bl = Bl4; Al = Al4; Bh = Bh4; Ah = Ah4;
      elseif(xframe_r_avg_dB <= -20)
      Wn = Xc(5);
      Bl = Bl5; Al = Al5; Bh = Bh5; Ah = Ah5;
    else
      Wn = Xc(6);
      Bl = Bl6; Al = Al6; Bh = Bh6; Ah = Ah6;
    end
    
    
    %% windowing
    
    xw_l = w .* xframe_l;                         % Apply window
    xw_r = w .* xframe_r;
    
    % HP for acoustic mic; LP for contact mic and filtering
    
    x_lp = filter(Bl, Al, xw_l);
    x_hp = filter(Bh, Ah, xw_r);
    
    % Combine LP and HP signals
    
    x_com = x_lp + x_hp;
    
end