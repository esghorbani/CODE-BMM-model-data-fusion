function [X_new] = TF(X, t_temp, dt)
  
    
    % Time integration using 4th-order Runge-Kutta
    
        Xt = X;
        
        k1 = BMM_UKF(t_temp, Xt);
        k2 = BMM_UKF(t_temp + dt/2, Xt + dt/2 * k1);
        k3 = BMM_UKF(t_temp + dt/2, Xt + dt/2 * k2);
        k4 = BMM_UKF(t_temp + dt, Xt + dt * k3);
        
        X_new = Xt + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    
end

