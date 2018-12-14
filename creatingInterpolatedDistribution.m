function outData=creatingInterpolatedDistribution(data,nBin,xq)

    [a1,b1]=hist(data,nBin);
    outData = interp1(b1,a1,xq);
        
end