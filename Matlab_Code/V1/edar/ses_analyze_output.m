function se_analyze_output(y,x,R,S)
% hFig = se_analyze_output(yout,Sout)

fs = R.fs ;

%hFig = figure ;
clf
nPlot = 4 ;
iPlot = 0 ;

%--------------------------------------------------------------------------
iPlot = iPlot+1 ;
ax(iPlot) = subplot(nPlot,1,iPlot) ;
timesignalplot(x(:,1),fs)
hold on
%timesignalplot(Sout.detect_result*max(abs(yout(:,1))),fs)
hold off
title('Contact microphone')

%--------------------------------------------------------------------------
iPlot = iPlot+1 ;
ax(iPlot) = subplot(nPlot,1,iPlot) ;
timesignalplot(x(:,2),fs)
hold on
%
hold off
title('Acoustic microphone')

%--------------------------------------------------------------------------
iPlot = iPlot+1 ;
ax(iPlot) = subplot(nPlot,1,iPlot) ;
timesignalplot(y(:,1),fs)
hold on
%
hold off
title(['Enhanced signal, WNS algo = ' S.WNS_method ', SE algo = ' S.SE_method])

%--------------------------------------------------------------------------
iPlot = iPlot+1 ;
ax(iPlot) = subplot(nPlot,1,iPlot) ;
timesignalplot(y(:,2),fs)
hold on
%
hold off
title('Reference signal (processed)')

%--------------------------------------------------------------------------
linkaxes(ax)

figure(gcf)



    function timesignalplot(Y,fs)
        hYline = plot(timevector(Y,fs),Y) ;
        grid on
        xlabel('Time (s)')
        %set(gca,'ButtonDownFcn',{@zoomplay,Y,fs})
%        [tb,btns] = axtoolbar({'zoomin','zoomout','restoreview'}) ;
%        btn = axtoolbarbtn(tb,'Tooltip','Play') ;
%        btn.ButtonPushedFcn = {@zoomplay,Y,fs} ;
        
    end

    function zoomplay(src,~,Y,fs)
        srcAxes = src.Parent.Parent ;
        axLim = axis(srcAxes) ;
        iLim = [max(round(axLim(1)*fs),1),...
                min(round(axLim(2)*fs),length(Y))] ;
        zoomGain = 1/max(abs(axLim(3:4))) ;
        if 0
            sound( zoomGain * Y(iLim(1):iLim(2)), fs)
        else
            player = audioplayer(zoomGain * Y(iLim(1):iLim(2)),fs) ;
            play(player)
            %disp('playing...')
            f = msgbox('Stop playback?','Audio Player','modal') ;
            while isplaying(player) && ishandle(f)
                pause(0.1)
            end
            stop(player)
            if ishandle(f)
                delete(f)
            end
            %disp('Stopped.')
        end
    end

end

