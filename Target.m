classdef Target <  handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        index
        position
        residingAgents
        
        % target state and error covariance 
        phi  % state 
        phiHat % state estimate
        Omega  % error covariuance of the state estimate
        
        timeSinceLastEvent;
        OmegaAtLastEvent;
        
        graphicHandles
        
        % From SAM
        A
        Q
        
        H
        R
        
        G
        
        % target control related paramters
        controllerEnabled
        B % control input matrix
        u % control input
        r % reference input
        K % controller gain
        
    end
    
    methods
        
        function obj = Target(index,A,Q,H,R,position)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.index = index;
            obj.position = position;
            
            % SAM
            obj.A = A;
            obj.Q = Q;
            
            obj.H = H;
            obj.R = R;
            
            obj.G = H'*inv(R)*H;
            

        end
        
        function outputArg = drawTarget(obj)
            hold on
            viscircles(obj.position,0.015,'Color','b');
            text(obj.position(1),obj.position(2)-0.04,num2str(obj.index),'Color','blue','FontSize',10);
        end
        
        function [cost,data,cost2] = update(obj,deltaT,plotMode,currentTime)
            
            obj.timeSinceLastEvent = obj.timeSinceLastEvent + deltaT;
            
            Omega1 = obj.Omega;
            controlError1 = abs(obj.phi-obj.r);
            
            % target state update
            w = mvnrnd(0,obj.Q); % w(t), zero mean with covariance Q, 
%                 obj.phi = obj.phi + deltaT*(obj.A*obj.phi + w); % phi_i(t+deltaT)
            obj.phi = obj.phi + deltaT*(obj.A*obj.phi + obj.B*obj.u + w); % phi_i(t+deltaT)
            
            % update of the state estimates
            if isempty(obj.residingAgents)
                obj.phiHat = obj.phiHat + deltaT*(obj.A*obj.phiHat + obj.B*obj.u); % phiHat_i(t+deltaT)
%                 obj.Omega = obj.Omega + deltaT*(obj.A*obj.Omega + obj.Omega*transpose(obj.A) + obj.Q); %Omega_i(t+deltaT)
                obj.updateLocalCovariance(0); % idle
            else
                v = mvnrnd(0,obj.R); % v(t), zero mean with covariance R, 
                z = obj.H*obj.phi + v; % z(t) observation 
                obj.phiHat = obj.phiHat + deltaT*(obj.A*obj.phiHat + obj.B*obj.u + obj.Omega*transpose(obj.H)*inv(obj.R)*(z-obj.H*obj.phiHat)); % phiHat_i(t+deltaT)
%                 obj.Omega = obj.Omega + deltaT*(obj.A*obj.Omega + obj.Omega*transpose(obj.A) + obj.Q - obj.Omega*obj.G*obj.Omega); %Omega_i(t+deltaT)
                obj.updateLocalCovariance(1); % dwell mode (active)
            end
            
            % control input computation (for the next time instant) based on the state estimate 
            if obj.controllerEnabled
                obj.r = 10*sin(2*currentTime+obj.index);
                rDot = 20*cos(2*currentTime+obj.index);
                obj.u = -(1/obj.B)*((obj.K+obj.A)*obj.phiHat - (rDot+obj.K*obj.r));
% %                 obj.u = -obj.K*obj.phiHat + (obj.K - obj.A/obj.B)*obj.r; % otherwise obj.u = 0 always
            else
                obj.u = 0;
            end
           
            % costs 
            Omega2 = obj.Omega;
            cost = (Omega1 + Omega2)*deltaT/2;
            cost2 = 0;
            if obj.controllerEnabled
                controlError2 = abs(obj.phi - obj.r);
                cost2 = (controlError1 + controlError2)*deltaT/2;
            end
            
            % data storage
            if plotMode
                eta_i = ~isempty(obj.residingAgents);
                if obj.controllerEnabled
                    data = [obj.phi, obj.phiHat, obj.Omega, obj.r, eta_i];
                else
                    data = [obj.phi, obj.phiHat, obj.Omega, 0, eta_i];
                end
            else
                data = [];
            end
            
        end
        
        function peakMin = computeLowerBound(obj,maxTime)
            maxTime = max(maxTime);
%             xss = (obj.A+sqrt(obj.A^2+obj.G*obj.Q))/obj.G;
%             peakMin = exp(2*obj.A*maxTime)*(xss + obj.Q*(1-exp(-2*obj.A*maxTime))/(2*obj.A));
            A_i = (obj.A + sqrt(obj.A^2 + obj.Q*obj.G))^2/(2*obj.A*obj.G);
            B_i = 2*obj.A;
            C_i = obj.Q/(2*obj.A);
            peakMin = A_i*exp(B_i*maxTime) - C_i;
        end
        
        
        function [gradientVal,intersectionVal] = computeGradientAndIntersection(obj, timeVal)
            
            A_i = (obj.A + sqrt(obj.A^2 + obj.Q*obj.G))^2/(2*obj.A*obj.G);
            B_i = 2*obj.A;
            C_i = obj.Q/(2*obj.A); 
            
            gradientVal = A_i*B_i*exp(B_i*timeVal);
            
            funVal = A_i*exp(B_i*timeVal) - C_i;
            intersectionVal = funVal - gradientVal*timeVal;
            
        end
        
        
        
        
        function J_i = localObjectiveIdlePerod(obj,startGap, duration) % startGap, duration;
            % omega_0 is the Omega after a dwell session of length startGap
            Omega_0 = obj.localCovarianceAfterDwellPerod(startGap);        
            a = obj.A;
            q = obj.Q;
            
            J_i = (1/(2*a))*(Omega_0 + (q/(2*a)))*(exp(2*a*duration)-1)-(q/(2*a))*duration;
        end
        
        
        function Omega = localCovarianceAfterDwellPerod(obj, duration)
            
            if duration==0
               Omega = obj.Omega;
               return;
            end
            
            a = obj.A;
            q = obj.Q;
            g = obj.G;
            v_1 = (-a+sqrt(a^2+q*g))/q;
            v_2 = (-a-sqrt(a^2+q*g))/q;
            lambda = 2*sqrt(a^2+q*g);
            Omega_0 = obj.Omega;
            
            c_1 = v_2*Omega_0 - 1;
            c_2 = -v_1*Omega_0 + 1;
            c_3 = v_1*c_1;
            c_4 = v_2*c_2;
            
            Omega = (c_1 + c_2*exp(-lambda*duration))/(c_3 + c_4*exp(-lambda*duration));
        end
        
        
        
        
        function J_i = localObjectiveDwellPerod(obj,startGap, duration) % startGap, duration;
            % omega_0 is the Omega after a idle session of length startGap
            Omega_0 = obj.localCovarianceAfterIdlePerod(startGap);
            a = obj.A;
            q = obj.Q;
            g = obj.G;
            v_1 = (-a+sqrt(a^2+q*g))/q;
            v_2 = (-a-sqrt(a^2+q*g))/q;
            lambda = 2*sqrt(a^2+q*g);
            
            c_1 = v_2*Omega_0 - 1;
            c_2 = -v_1*Omega_0 + 1;
            c_3 = v_1*c_1;
            c_4 = v_2*c_2;
            
%             J_i = (1/g)*log(abs(c_3+c_4*exp(-lambda*duration))) + (1/v_1)*duration - (1/g)*log(abs(v_2-v_1));
            J_i = (1/g)*log((c_3+c_4*exp(-lambda*duration))/(v_2-v_1)) + (1/v_1)*duration;
        end
        
        
        function Omega = localCovarianceAfterIdlePerod(obj, duration)
            
            if duration==0
               Omega = obj.Omega;
               return;
            end
            
            a = obj.A;
            q = obj.Q;
            Omega_0 = obj.Omega;
            
            Omega = (Omega_0 + (q/(2*a)))*exp(2*a*duration) - q/(2*a);
        end
        
        
        function dwellTime = dwellTimeToReduceCovarianceUptoFraction(obj, e)
        
            a = obj.A;
            q = obj.Q;
            g = obj.G;
            v_1 = (-a+sqrt(a^2+q*g))/q;
            v_2 = (-a-sqrt(a^2+q*g))/q;
            lambda = 2*sqrt(a^2+q*g);
            
            Omega_0 = obj.Omega;
            c_1 = v_2*Omega_0 - 1;
            c_2 = -v_1*Omega_0 + 1
%             c_3 = v_1*c_1;
%             c_4 = v_2*c_2;
            
            Omega_ss = 1/v_1;
            if Omega_ss < Omega_0 - 0.001 & c_2 ~= 0
                dwellTime = (-1/lambda)*log((-e*v_1*c_1)/((v_1+(1-e)*v_2)*c_2));
                if dwellTime>-2
                    dwellTime = abs(dwellTime);
                end
            elseif Omega_ss > Omega_0
%                 dwellTime = (-1/lambda)*log((-e*v_1*c_1)/((v_1+(1-e)*v_2)*c_2));
                dwellTime = 0;
            else
                dwellTime = 0;
            end
            
            if dwellTime<0
                disp(['Target ',num2str(obj.index),' observation already complete.']);
                dwellTime = 0;
            elseif imag(dwellTime)~=0
                disp(['Target ',num2str(obj.index),' observation error.']);
                dwellTime = 0;
            end
                
            
        
        end
        
        function output = updateLocalCovariance(obj, mode)
            if mode == 0 % idle
                a = obj.A;
                q = obj.Q;
                Omega_0 = obj.OmegaAtLastEvent;
                duration = obj.timeSinceLastEvent;
                obj.Omega = (Omega_0 + (q/(2*a)))*exp(2*a*duration) - q/(2*a);
            else % dwell
                a = obj.A;
                q = obj.Q;
                g = obj.G;
                v_1 = (-a+sqrt(a^2+q*g))/q;
                v_2 = (-a-sqrt(a^2+q*g))/q;
                lambda = 2*sqrt(a^2+q*g);
                Omega_0 = obj.OmegaAtLastEvent;
                duration = obj.timeSinceLastEvent;

                c_1 = v_2*Omega_0 - 1;
                c_2 = -v_1*Omega_0 + 1;
                c_3 = v_1*c_1;
                c_4 = v_2*c_2;
                
                obj.Omega = (c_1 + c_2*exp(-lambda*duration))/(c_3 + c_4*exp(-lambda*duration));
            end
        end
        
        
    end
end

