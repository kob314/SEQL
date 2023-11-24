function [ks kl] = ks_kl_fun(tt,err)

    kl_aux = 1./(pi*sqrt(tt));
    kl = max( sqrt(-2*log(pi*tt*err)./(pi^2*tt)) , kl_aux ); % bound
    
    kl((pi*tt*err)>=1) = kl_aux((pi*tt*err)>=1);

    ks = max( 2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)) , sqrt(tt)+1 );
    ks((2*sqrt(2*pi*tt)*err)>=1)=2;

