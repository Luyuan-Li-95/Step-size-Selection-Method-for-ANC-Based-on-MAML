%% Modified MAML algorithm
classdef MAML_train_test_forget
    properties 
        Phi % The initial step size 
        W_star %one step optimal control filter
    end
    methods 
        function obj = MAML_train_test_forget(initial,len_c)
            % give stepsize an initial value
            obj.Phi = initial;
            obj.W_star = zeros(len_c,1);
        end
        function [obj,Er] = MAML_initial(obj,Fx,Di,lamda,epslon)
            % Fx : the filtered reference vector 
            % Di : the disturbance vector  
            % W  : the initial value of the control filter 
            % lamda : the forget factor 
            Fx   = flipud(Fx);
            Dis  = flipud(Di); 
            Grad = 0; % Temporal gradient accumulator
            Er   = 0; % Training error signal 
            Li   = length(obj.W_star) ; % The length of the control filter in the FxLMS algorithm. 
			%<-4-> Get the error signal. 
			e    = Dis(1)  - obj.W_star'*Fx; 
			%<-5-> Obtain the control filter
            Wo   = obj.W_star + obj.Phi*e*Fx    ; % One-step updation for the assumed optimal control filter.
            for jj = 1:Li
                Fd   = [Fx(jj:end);zeros(jj-1,1)];
				%<-6-> Get the error signal based on the new control filter.
                e_dagger    = Dis(jj) - Wo'*Fd          ; 
				% Get the gradints based on the assumed optimal control filter 
                % Grad = Grad    + epslon*(mu/Li)*e*Fd*(lamda^(jj-1)); 
                Grad = Grad    + epslon/Li*(lamda^(jj-1))*e_dagger*e*Fd'*Fx; 
                if jj == 1
                    Er =  e_dagger ;
                end
            end
			%%<-7-> Upate the initial value 
            obj.Phi = obj.Phi + Grad ;
            obj.W_star = Wo;
        end
    end
    
end
%------------------------- end --------------------------------------------