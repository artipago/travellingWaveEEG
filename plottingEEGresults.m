function plottingEEGresults(dataEEG,dataEEGNull)


m=20;
xq=-1.7:0.05:1.7;
logRatioInput=[];
logRatioInputNull=[];
fwValue=[];
bwValue=[];
fwValueSS=[];
bwValueSS=[];
binsHist=3:20;
for ss=1:size(dataEEG,2) %here we combine all the sbj in one vector
    logRatioInput=[logRatioInput dataEEG(ss).logRatio];
    logRatioInputNull=[logRatioInputNull dataEEGNull(ss).logRatio];
    
    fwValue=[fwValue nanmean(dataEEG(ss).fwValue)];
    bwValue=[bwValue nanmean(dataEEG(ss).bwValue)];
    
    fwValueSS=[fwValueSS nanmean(dataEEGNull(ss).fwValue)];
    bwValueSS=[bwValueSS nanmean(dataEEGNull(ss).bwValue)];
    
    fwTempFreq(ss,:)=gettingTheHist(dataEEG(ss).fwTempFreq,binsHist);
    bwTempFreq(ss,:)=gettingTheHist(dataEEG(ss).bwTempFreq,binsHist);
    fwTempFreqSS(ss,:)=gettingTheHist(dataEEGNull(ss).fwTempFreq,binsHist);
    bwTempFreqSS(ss,:)=gettingTheHist(dataEEGNull(ss).bwTempFreq,binsHist);
end

inputDist=creatingInterpolatedDistribution(logRatioInput,m,xq);
inputDistNull=creatingInterpolatedDistribution(logRatioInputNull,m,xq);

[input,inputSS,diffInput,inputFW,inputBW]=normalizeNcomputeDiff(inputDist,inputDistNull,xq);

figure
subplot(2,2,1)
hold on
plot(xq,input,'k')
plot(xq,inputSS,'--k')
title('Distribution of logRatios')
legend('real distribution','null distribution')
xlabel('logRatio')
xlim([min(xq) max(xq)])

subplot(2,2,3)
hold on
plot(xq,diffInput,'k')
title('Difference between real and null distributions')
xlabel('logRatio')
xlim([min(xq) max(xq)])

subplot(3,2,2)
hold on
errorbar(1.1,mean(bwValue),std(bwValue)./sqrt(length(bwValue)),'r')
errorbar(1.9,mean(fwValue),std(fwValue)./sqrt(length(fwValue)),'b')
errorbar(1.4,mean(bwValueSS),std(bwValueSS)./sqrt(length(bwValueSS)),'k')
errorbar(1.6,mean(fwValueSS),std(fwValueSS)./sqrt(length(fwValueSS)),'k')
plot(1,bwValue,'r*')
plot(2,fwValue,'b*')
xlim([0.5 2.5])
legend('BW','FW','Null')
title('Max power value from 2DFFT')

subplot(3,2,4)
hold on
plot(binsHist,mean(bwTempFreq),'r')
plot(binsHist,mean(fwTempFreq),'b')
plot(binsHist,mean(bwTempFreqSS),'k')
plot(binsHist,mean(fwTempFreqSS),'k')
% plot(1,bwTempFreq,'r*')
% plot(2,fwTempFreq,'b*')
legend('BW','FW','Null')
title('Temporal frequency')

subplot(3,2,6)
hold on
bar(1,inputFW,'b')
bar(2,inputBW,'r')
title('percentage of significant wave')
legend('FW','BW')
xlim([0.5 2.5])



