function [logRatio,bwValue,bwTempFreq,bwSpatFreq,fwValue,fwTempFreq,fwSpatFreq]=wavesHunter(data,samplingRate,freqBandFlag,simmetryFlag)

%% data has size: [electrode time]. freqBandFlag and simmetryFlag defines how we analyze the 2DFFT. 
% if freqBandFlag==1 it searches in the quadrant from 2Hz to 100Hz and pick the maximum value. 
% if freqBandFlag==2 it searches in the quadrant from 8Hz to 12Hz and pick the maximum value there (focus in the alpha-band). 
% if simmetryFlag==1 it takes the maximum value over both quadrant, and its symmetric in the other quadrant. 
% if simmetryFlag==0 it takes maximum values in each quadrant, indepdently of the other (no simmetry). 

    numberElectrodes=size(data,1);
    yElec=1:numberElectrodes;
    durationSignal=size(data,2)/samplingRate; %it's in second.
    dF = 1/durationSignal;                   
    fx = -samplingRate/2:dF:samplingRate/2-dF;         
    fy=(size(data,1)/2)*linspace(-1,1,size(data,1));
    sbj2DFFT=abs(fftshift(fft2(data)));
    
    if freqBandFlag==1
        f1=2;
        f2=100;
    elseif freqBandFlag==2
        f1=7;
        f2=13;
    else
        error('freqBandFlag should be either 1 (all spectra) or 2 (only alphaBand)')
    end
    
    if simmetryFlag==0
        [aBW,bBW]=max(sbj2DFFT(fy>0, fx>=f1 & fx<=f2));
        [bwValue,bBW2]=max(aBW);
        bwTempFreq=fx(bBW2+ceil(length(fx)/2)+round(f1/dF));
        bwSpatFreq=fy(bBW(bBW2)+ceil(length(yElec)/2));

        
        [aFW,bFW]=max(sbj2DFFT(fy<0, fx>=f1 & fx<=f2));
        [fwValue,bFW2]=max(aFW);
        fwTempFreq=fx(bFW2+ceil(length(fx)/2)+round(f1/dF));
        fwSpatFreq=fy(bFW(bFW2));
    else
        if mod(length(yElec),2)==1 %if odd it has a zero spatial component to kill when looking for the general maximum
            yElec(ceil(length(yElec)/2))=[]; 
        end
        [a,b]=max(sbj2DFFT(yElec, fx>=f1 & fx<=f2));
        b(b>3)=b(b>3)+1; %to keep into account the zeroSpatial component
        [a2,b2]=max(a);
        if b(b2)<=ceil(length(yElec)/2) 
            %maximuml is in the FW
            fwValue=a2;
            bwValue= sbj2DFFT(-b(b2)+numberElectrodes+1,b2+ceil(length(fx)/2)+round(f1/dF));
            fwSpatFreq=fy(b(b2));
            bwSpatFreq=fy(-b(b2)+numberElectrodes+1);
        else
            %maximuml is in the BW
            bwValue=a2;
            fwValue= sbj2DFFT(-b(b2)+numberElectrodes+1,b2+ceil(length(fx)/2)+round(f1/dF));

            bwSpatFreq=fy(b(b2));
            fwSpatFreq=fy(b(b2)-(ceil(length(yElec)/2))); %or +1
        end
        bwSpatFreq=abs(fy(b(b2)));
        fwSpatFreq=-abs(fy(b(b2)));
        bwTempFreq=fx(b2+ceil(length(fx)/2)+round(f1/dF));
        fwTempFreq=fx(b2+ceil(length(fx)/2)+round(f1/dF));
    end
    
    logRatio=log(fwValue/bwValue);


    

end
