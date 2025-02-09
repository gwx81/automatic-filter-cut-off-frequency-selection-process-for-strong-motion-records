function [tf]=fbutter(Nnyq,df,fc,n,hp,AC)
if AC=='A'
    f=df*((1:Nnyq)-1);
    if hp=='h'
        tf=sqrt((f./fc).^(2*n)./(1.0+(f./fc).^(2*n)));
    end
    if hp=='l'
        tf = sqrt(1.0./(1.0+(f./fc).^(2*n)));
    end
end

if AC=='C'
    f=df*((1:Nnyq)-1);
    tf = complex(1,0);
    if hp == 'h'
        fn = -complex(0,abs(fc./f));
        for j = 1:n
            tf = tf.*(fn-exp(complex(0,pi/(2.*n)*(2.*j-1+n))));
        end
        tf = 1./tf;
    end
    if hp == 'l'
        fn = complex(0,abs(f./fc));
        for j= 1:n
            tf = tf .* (fn-exp(complex(0,pi/(2.*n)*(2.*j-1+n))));
        end
        tf = 1./tf;
    end
end
end