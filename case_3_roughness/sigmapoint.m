function [W,SP] = sigmapoint(X,p,gama,L)
% Unscented Kalman Filter Code 

     W(1)=gama/(L+gama);  % sigma point weight
     
     eps_reg = 1e-19 * eye(size(p));
     d=chol((L+gama)*(p + eps_reg));  % squre root of P matrix caclulation 
     
     SP(:,1)=X;    % Initial Sigma Point
     
     % Sigma points and its related weight value calculation
     
     for c=1:L
     SP(:,c+1)=SP(:,1)+d(c,:)';
     W(c+1)=1/(2*(L+gama));
     end
     
         
     for j=1:L
     SP(:,j+L+1)=SP(:,1)-d(j,:)';
     W(j+L+1)=1/(2*(L+gama));

end

