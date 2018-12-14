function [f1new,f2new,diffF,fw,bw]=normalizeNcomputeDiff(f1,f2,xq)

    f1new=f1./nansum(f1);
    f2new=f2./nansum(f2);
    diffF=f1new-f2new;
    
    normDiff=diffF./nansum(abs(diffF));
    normDiff(normDiff<0)=0;
    fw=nansum(normDiff(xq>0));
    bw=nansum(normDiff(xq<0));
    
    diffF(diffF<0)=0;
    
end
