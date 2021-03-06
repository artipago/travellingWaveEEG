function [vOut,fw,bw]=binningVect(vIn,nBins)

xx=1:floor(length(vIn)/nBins):length(vIn);
for ii=1:nBins-1
    tempV=vIn(xx(ii):xx(ii+1));
    tempV(tempV==0)=NaN;
    vOut(ii)=nanmean(tempV);
    fw(ii)=sum(vIn(xx(ii):xx(ii+1))>0);
    bw(ii)=sum(vIn(xx(ii):xx(ii+1))<0);
end
vOut(nBins)=mean(vIn(xx(nBins:end)));
% fw(nBins)=sum(vIn(xx(nBins:end))>0);
% bw(nBins)=sum(vIn(xx(nBins:end))<0);
