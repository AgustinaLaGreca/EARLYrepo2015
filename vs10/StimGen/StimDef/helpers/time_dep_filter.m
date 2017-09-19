function fsignal = time_dep_filter(signal,fcs,width,order,Fs,type,interval)

nsamples = length(signal);
fsignal = zeros(1,nsamples);
make_coeff=1;


for ifc=1:nsamples-order
%     ifc,nsamples-order
    fc=fcs(ifc);
    
    f1 = fc/(2^(width/2));
    f2= fc*(2^(width/2));
    Wn = [f1/Fs*2,f2/Fs*2];
%     fc,f1,f2
    
    if ifc==1
        b = fir1(order,Wn,type);
        b=b/sum(b);
    end
    if make_coeff==interval
%         ifc,interval
        b = fir1(order,Wn,type);
        b=b/sum(b);
        make_coeff=1;
    else
        make_coeff=make_coeff+1;
    end
    fsignal(ifc)=signal(ifc:ifc+order)*b(end:-1:1)';
end
end