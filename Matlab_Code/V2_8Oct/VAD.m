
function [vad_signal, vad_average] = VAD(xframe_l_filtered, vad_average, alpha, M, nhop, vad_threshold );
    p = vad_average;
    for k=1:M
      p = (1-alpha)*p+alpha*xframe_l_filtered(k).^2;
      vad_average(k) = p;
      
      if (E >= vad_threshold)
        vad_signal(k) = 1;
      else
        vad_signal(k) = 0;
      endif
    endfor
endfunction
