function X_new = TF_rough(X, t_temp, dt, par)
%TF_ROUGH  One RK4 step of the estimator model BMM_UKF_ROUGH (cf. TF.m).

k1 = BMM_UKF_rough(t_temp,        X,             par);
k2 = BMM_UKF_rough(t_temp + dt/2, X + dt/2 * k1, par);
k3 = BMM_UKF_rough(t_temp + dt/2, X + dt/2 * k2, par);
k4 = BMM_UKF_rough(t_temp + dt,   X + dt   * k3, par);

X_new = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
