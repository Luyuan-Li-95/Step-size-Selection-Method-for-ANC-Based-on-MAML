%% single-channel FxNLMS algorithm 
function [Er, muw] = FxNLMS(Len_Filter, Wc_initial, Dis, Rf, Beta)
% Len_Filter : the length of the control filter 
% Wc_initial : the initial control filter 
% Dis        : the disturbance 
% Rf         : the filtered reference vector 
% Beta       : the parameter for adjusting the step size 
N   = Len_Filter ;
Wc  = Wc_initial ;
XD  = zeros(N,1) ;
Er  = zeros(length(Rf),1);
Pe  = zeros(length(Rf),1);
muw = zeros(length(Rf),1);
    for tt = 1:length(Rf) 
        XD   = [Rf(tt);XD(1:end-1)];
        Rf_i = XD'         ;
        Rf_i = Rf_i'       ;
        y_t  = Wc'*Rf_i    ;
        e    = Dis(tt)-y_t ;
        Er(tt) = e         ;
        if tt==1
            Pe(tt) = (norm(XD))^2;
        else
            Pe(tt) = 0.9*Pe(tt-1) + 0.1*(norm(XD))^2;
        end
        muw(tt)    = Beta/Pe(tt);
        if muw(tt)>1e-3
            muw(tt) = 1e-3;
        else
            muw(tt) = muw(tt);
        end

        Wc     = Wc + muw(tt)*e*Rf_i;
    end
end
%-------------------------end-----------------------------