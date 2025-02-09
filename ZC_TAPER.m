function [ACC]=ZC_TAPER(acc,DT)
TB=0;
TE=0;
for ii=2:length(acc)
    if acc(ii,1)*acc(ii-1,1)<0 && ii>=5
        TB=ii*DT;
        break
    end
end
if ii==length(acc)
    TB=0;
    acc(1,1)=0;
end
for ii=length(acc):-1:2
    if acc(ii,1)*acc(ii-1,1)<0 && ii<=length(acc)-0.01*length(acc)
        TE=(length(acc)-ii)*DT;
        break
    end
end
if ii==2
    TE=0;
    acc(length(acc),1)=0;
end
ACC(:,1)=Taper(DT,acc,TB,TE);

end