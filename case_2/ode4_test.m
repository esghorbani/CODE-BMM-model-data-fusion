function [t, X,Mt,Kt,Ct,ft] = ode4_test(odefun, tspan, X0, dt)
    % Number of steps
    nSteps = floor((tspan(end) - tspan(1)) / dt);
    
    % Initialize the solution
    t = linspace(tspan(1), tspan(end), nSteps+1)';
    X = zeros(nSteps+1, length(X0));
    X(1,:) = X0';
    
    % Time integration using 4th-order Runge-Kutta
    for i = 1:nSteps
        tt = t(i);
        Xt = X(i,:)';
        
        k1 = odefun(tt, Xt);
        k2 = odefun(tt + dt/2, Xt + dt/2 * k1);
        k3 = odefun(tt + dt/2, Xt + dt/2 * k2);
        [k4,Mt,Kt,Ct,ft] = odefun(tt + dt, Xt + dt * k3);
        
        X(i+1,:) = Xt + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end

