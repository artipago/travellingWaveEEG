function plottingEchoResults(dataEcho,dataEchoNull)

m=20;
xq=-1.7:0.05:1.7;

echoesDist=creatingInterpolatedDistribution(dataEcho.logRatio,m,xq);
echoesDistNull=creatingInterpolatedDistribution(dataEchoNull.logRatio,m,xq);

[echoes,echoesSS,diffEcho,echoesFW,echoesBW]=normalizeNcomputeDiff(echoesDist,echoesDistNull,xq);

figure
subplot(2,2,1)
hold on
plot(xq,echoes,'k')
plot(xq,echoesSS,'--k')
title('Distribution of logRatios')
legend('real distribution','null distribution')
xlabel('logRatio')
xlim([min(xq) max(xq)])

subplot(2,2,3)
hold on
plot(xq,diffEcho,'k')
title('Difference between real and null distributions')
xlabel('logRatio')
xlim([min(xq) max(xq)])

subplot(3,2,2)
hold on
errorbar(1.1,mean(dataEcho.bwValue),std(dataEcho.bwValue)./sqrt(length(dataEcho.bwValue)),'r')
errorbar(1.9,mean(dataEcho.fwValue),std(dataEcho.fwValue)./sqrt(length(dataEcho.fwValue)),'b')
errorbar(1.4,mean(dataEchoNull.bwValue),std(dataEchoNull.bwValue)./sqrt(length(dataEchoNull.bwValue)),'k')
errorbar(1.6,mean(dataEchoNull.fwValue),std(dataEchoNull.fwValue)./sqrt(length(dataEchoNull.fwValue)),'k')
plot(1,dataEcho.bwValue,'r*')
plot(2,dataEcho.fwValue,'b*')
xlim([0.5 2.5])
axis([0.5 2.5 0 5])
legend('BW','FW','Null')
title('Max power value from 2DFFT')

subplot(3,2,4)
hold on
errorbar(1.2,mean(dataEcho.bwTempFreq),std(dataEcho.bwTempFreq)./sqrt(length(dataEcho.bwTempFreq)),'r')
errorbar(1.8,mean(dataEcho.fwTempFreq),std(dataEcho.fwTempFreq)./sqrt(length(dataEcho.fwTempFreq)),'b')
errorbar(1.4,mean(dataEchoNull.bwTempFreq),std(dataEchoNull.bwTempFreq)./sqrt(length(dataEchoNull.bwTempFreq)),'k')
errorbar(1.6,mean(dataEchoNull.fwTempFreq),std(dataEchoNull.fwTempFreq)./sqrt(length(dataEchoNull.fwTempFreq)),'k')
% plot(1,dataEcho.bwTempFreq,'r*')
% plot(2,dataEcho.fwTempFreq,'b*')
axis([1 2 5 12])
legend('BW','FW','Null')
title('Temporal frequency')

subplot(3,2,6)
hold on
bar(1,echoesFW,'b')
bar(2,echoesBW,'r')
title('percentage of significant wave')
legend('FW','BW')
xlim([0.5 2.5])






