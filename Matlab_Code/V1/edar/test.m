clc
close all
clear all
Wn = 5000/(24000) ;             % Normalized crossover frequency
[Bl,Al] = butter(2,Wn,'low')      % Butterworth filter
[Bh,Ah] = butter(2,Wn,'high')     % Butterworth filter



freqz(Bl,Al);