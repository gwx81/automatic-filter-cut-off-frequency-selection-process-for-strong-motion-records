function[acctaper]=Taper(dt,acc,tbegin,tend)


if tbegin~=0 && tend ==0
    nb=floor(tbegin/dt);
    for i=1:nb
        resultb(i,1)=1/2*(1+cos(pi*(nb+i-1)/nb));
    end
    ne=1;
    resulte=1;
    
elseif tend~=0 && tbegin ==0
    ne=floor(tend/dt);
    for i=1:ne
        resulte(i,1)=1/2*(1+cos(pi*(i-1)/ne));
    end
    nb=1;
    resultb=1;
    
elseif tend~=0 && tbegin ~=0
    ne=floor(tend/dt);
    for i=1:ne
        resulte(i,1)=1/2*(1+cos(pi*(i-1)/ne));
    end
    
    nb=floor(tbegin/dt);
    for i=1:nb
        resultb(i,1)=1/2*(1+cos(pi*(nb+i-1)/nb));
    end
else
    nb=1;
    resultb=1;
    ne=1;
    resulte=1;
    
end

acc1(:,1)=acc(1:nb,1).*resultb(:,1);
acc2(:,1)=acc(end-ne+1:end,1).*resulte(:,1);
acctaper(:,1)=[acc1(:,1)' acc(nb+1:end-ne,1)' acc2(:,1)']';