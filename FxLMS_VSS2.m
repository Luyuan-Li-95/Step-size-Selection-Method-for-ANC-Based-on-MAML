%%  
function [Er, muw] = FxLMS_VSS2(Len_Filter, Wc_initial, Dis, Rf, Beta)
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
Px  = zeros(length(Rf),1);
muw = zeros(length(Rf),1);
    for tt = 1:length(Rf) 
        XD   = [Rf(tt);XD(1:end-1)];
        Rf_i = XD'         ;
        Rf_i = Rf_i'       ;
        y_t  = Wc'*Rf_i    ;
        e    = Dis(tt)-y_t ;
        Er(tt) = e         ;
        if tt==1
            Px(tt) = 0.05*(norm(XD))^2;
            Pe(tt) = 0.05*e^2;
        else
            Px(tt) = 0.95*Px(tt-1) + 0.05*(norm(XD))^2;
            Pe(tt) = 0.95*Pe(tt-1) + 0.05*e^2;
        end
        muw(tt)    = Beta / (sqrt(Px(tt))+Pe(tt)+1e-8);
        if muw(tt)>1e-3
            muw(tt) = 1e-3;
        else
            muw(tt) = muw(tt);
        end

        Wc     = Wc + muw(tt)*e*Rf_i;
    end
end
%-------------------------end-----------------------------