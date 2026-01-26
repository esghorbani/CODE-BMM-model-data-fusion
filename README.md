# BMM_model_data_fusion

This repository contains the reference implementation for our paper (submitted to Mechanical Systems and Signal Processing, MSSP) on joint mass–velocity estimation of a moving object traversing an Euler–Bernoulli beam from single-point vibration measurements using an Adaptive Unscented Kalman Filter (AUKF).

We formulate the moving object as a traveling mass that induces inertial, centrifugal, and Coriolis interaction forces, and embed the governing equations into an augmented state-space model where mass and velocity are treated as unknown parameters. The UKF is derived from Kalman filtering principles, and the measurement noise covariance is tuned online via a residual-based adaptive scheme. In addition to parameter identification, the framework reconstructs the full dynamic contact force history, enabling applications such as bridge weight-in-motion (WIM) without requiring an external speed sensor.

The code reproduces the main experiments reported in the paper, including the comparison of displacement vs. acceleration measurements for robustness and convergence across a wide range of mass–velocity combinations.
