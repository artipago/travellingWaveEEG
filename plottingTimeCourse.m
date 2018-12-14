function plottingTimeCourse(dataEEG,dataEEGnull)

nBins=60;
for ss=1:size(dataEEG,2)
    q025=quantile(dataEEGnull(ss).logRatio,0.025);
    q975=quantile(dataEEGnull(ss).logRatio,0.975);
    logRatio=dataEEG(ss).logRatio;
    logRatio(logRatio>q025 & logRatio<q975)=0;
%     quantiles(ss,:)=[q025 q975];
    [timeCourse(ss,:),fw(ss,:),bw(ss,:)]=binningVect(logRatio,nBins);
end

figure
subplot(2,1,1)
hold on

errorbar(nanmean(timeCourse,1),nanstd(timeCourse,0,1)./sqrt(size(timeCourse,1)))
plot([1 nBins],[0 0],'k')
title('trends of LogRatios above&below 0.05 quantile of the null distribution')
ylabel('logRatios')
xlabel('time - in bins')

subplot(2,1,2)
hold on
plot(sum(fw,1),'b')
plot(sum(bw,1),'r')
ylabel('waves count')
legend('BW','FW')
xlabel('time - in bins')

figure
for ss=1:size(timeCourse,1)
    subplot(ceil(size(timeCourse,1)/6),6,ss)
    hold on
    plot([0 nBins],[0 0],'r')
    plot(timeCourse(ss,:),'b')
%     axis([0 nBins -1.5 1.5])
    title(['sbj ' int2str(ss)])
    ylabel('logRatios')
    xlabel('time ')
end
