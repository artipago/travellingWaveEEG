function res=gettingTheHist(vect,bins)

    for ii=2:length(bins)
        res(ii)=length(vect(vect>bins(ii-1) & vect<=bins(ii)));
    end
    
end