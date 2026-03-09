%%  
function [Er, muw] = FxLMS_CSS(Len_Filter, Wc_initial, Dis, Rf, beta, mu1, mu2)
% Len_Filter : the length of the control filter 
% Wc_initial : the initial control filter 
% Dis        : the disturbance 
% Rf         : the filtered reference vector 
% beta       : the parameter for adjusting the step size 
% mu1,mu2    : 2 step size
lamda = 0.99;
zeta = 0.9;
mub = 0.1;

N   = Len_Filter ;
Wc  = Wc_initial ;
XD  = zeros(N,1) ;
Er  = zeros(length(Rf),1);
g   = zeros(length(Rf),1);
b   = zeros(length(Rf),1);
k   = zeros(N,length(Rf));
gama= zeros(length(Rf),1);
muw = zeros(length(Rf),1);
    for tt = 1:length(Rf) 
        XD   = [Rf(tt);XD(1:end-1)];
        Rf_i = XD'         ;
        Rf_i = Rf_i'       ;
        y_t  = Wc'*Rf_i    ;
        e    = Dis(tt)-y_t ;
        Er(tt) = e         ;
        
        % calculate step size
        if tt==1
            g(tt) = Dis(tt);
            phi = eye(N,N) .* (0.04)^(-1);
            b(tt) = 0;
            gama(tt) = 1;
        else
            g(tt) = Dis(tt) - XD'*Wc_d;
            m     = 1 / (1+beta*(g(tt)^2))^2;
            u     = phi * XD;
            k(:,tt) = u / (lamda/m + XD'*u);
            phi = zeta.*phi + (1-zeta)*lamda^(-1).*(phi-k(:,tt)*u');
            b(tt) = b(tt-1) + mub * sign(g(tt-1)) * g(tt-1) * gama(tt-1)*(1-gama(tt-1)) .* XD'*k(:,tt-1);
            if b(tt)<-2.49
                b(tt) = -2.49;
            elseif b(tt)>2.49
                b(tt) = 2.49;
            else
                b(tt) = b(tt);
            end
            gama(tt) = 0.5 / (1+exp(-b(tt)));
        end
        muw(tt) = gama(tt)*mu1 + (1-gama(tt))*mu2;
        
        Wc_d   = Wc;
        Wc     = Wc + muw(tt)*e*Rf_i;
    end
end
%-------------------------end-----------------------------