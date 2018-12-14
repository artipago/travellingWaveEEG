function startingFunction()

%% loading example data

% echoesInput, eegInput, eegClosed. All are structures of sbj (n=20 and
% 48). Sampling rate = 160Hz.

% Data sizes are: 
% echoesInput(ss).data = [electrode time] 
% eegInput(ss).data = [electrode time trials] 
% eegClosed(ss).data = [electrode time] 

load('dataExample.mat') 


%% analysis echoes

samplingRate=160;
freqBandFlag=1; %if 1 consider all spectra from 2Hz, if 2 only within the alphaBand
simmetryFlag=0; %if 1 the maxima are simmetric (the lower may not be the maximum of its quadrant), if 0 they are independent.

for ss=1:size(echoesInput,2)
    [dataEcho.logRatio(ss),dataEcho.bwValue(ss),dataEcho.bwTempFreq(ss),dataEcho.bwSpatFreq(ss),dataEcho.fwValue(ss),dataEcho.fwTempFreq(ss),dataEcho.fwSpatFreq(ss)]=wavesHunter(echoesInput(ss).data,samplingRate,freqBandFlag,simmetryFlag);
    [dataEchoNull.logRatio(ss),dataEchoNull.bwValue(ss),dataEchoNull.bwTempFreq(ss),dataEchoNull.bwSpatFreq(ss),dataEchoNull.fwValue(ss),dataEchoNull.fwTempFreq(ss),dataEchoNull.fwSpatFreq(ss)]=wavesHunter(echoesInput(ss).data(randperm(size(echoesInput(ss).data,1)),:),samplingRate,freqBandFlag,simmetryFlag);
end

plottingEchoResults(dataEcho,dataEchoNull)


%% analysis eegInput

samplingRate=160;
freqBandFlag=1; %if 1 consider all spectra from 2Hz, if 2 only within the alphaBand
simmetryFlag=1; %if 1 the maxima are simmetric (the lower may not be the maximum of its quadrant), if 0 they are independent.

for ss=1:size(eegInput,2) %running over sbj
    for tt=1:size(eegInput(ss).data,3) %running over trials
        [dataInput(ss).logRatio(tt),dataInput(ss).bwValue(tt),dataInput(ss).bwTempFreq(tt),dataInput(ss).bwSpatFreq(tt),dataInput(ss).fwValue(tt),dataInput(ss).fwTempFreq(tt),dataInput(ss).fwSpatFreq(tt)]=wavesHunter(eegInput(ss).data(:,:,tt),samplingRate,freqBandFlag,simmetryFlag);
        [dataInputNull(ss).logRatio(tt),dataInputNull(ss).bwValue(tt),dataInputNull(ss).bwTempFreq(tt),dataInputNull(ss).bwSpatFreq(tt),dataInputNull(ss).fwValue(tt),dataInputNull(ss).fwTempFreq(tt),dataInputNull(ss).fwSpatFreq(tt)]=wavesHunter(eegInput(ss).data(randperm(size(eegInput(ss).data,1)),:,tt),samplingRate,freqBandFlag,simmetryFlag);
    end
end

plottingEEGresults(dataInput,dataInputNull)
plottingTimeCourse(dataInput,dataInputNull)


%% analysis eegClosed
samplingRate=1024;
sizeTimeWindow=1; % in seconds
overlapBetweenTW=0.5; %in seconds
freqBandFlag=1; %if 1 consider all spectra from 2Hz, if 2 only within the alphaBand
simmetryFlag=0; %if 1 the maxima are simmetric (the lower may not be the maximum of its quadrant), if 0 they are independent.

for ss=1:size(eegClosed,2) %running over sbj
    xx=1:overlapBetweenTW*samplingRate:size(eegClosed(ss).data,2)-overlapBetweenTW*samplingRate; %starting point of each time window
    for ii=1:length(xx)-1 %running over TW
        [dataClosed(ss).logRatio(ii),dataClosed(ss).bwValue(ii),dataClosed(ss).bwTempFreq(ii),dataClosed(ss).bwSpatFreq(ii),dataClosed(ss).fwValue(ii),dataClosed(ss).fwTempFreq(ii),dataClosed(ss).fwSpatFreq(ii)]=wavesHunter(eegClosed(ss).data(:,xx(ii):xx(ii)+sizeTimeWindow*samplingRate),samplingRate,freqBandFlag,simmetryFlag);
        [dataClosedNull(ss).logRatio(ii),dataClosedNull(ss).bwValue(ii),dataClosedNull(ss).bwTempFreq(ii),dataClosedNull(ss).bwSpatFreq(ii),dataClosedNull(ss).fwValue(ii),dataClosedNull(ss).fwTempFreq(ii),dataClosedNull(ss).fwSpatFreq(ii)]=wavesHunter(eegClosed(ss).data(randperm(size(eegClosed(ss).data,1)),xx(ii):xx(ii)+sizeTimeWindow*samplingRate),samplingRate,freqBandFlag,simmetryFlag);
    end
end

plottingEEGresults(dataClosed,dataClosedNull)
plottingTimeCourse(dataClosed,dataClosedNull)


