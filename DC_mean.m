function [acc1]=DC_mean(acc)
    average=mean(acc);
    acc1=acc-average;
end