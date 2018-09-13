kalmanBased-eventDetection
A threshold-free, preprocessing-free, fast Kalman filter-based gait event detection system. 
It takes forces along the gravity axis for each individual leg, and estimates the stance/swing phases for each leg (and consequently the relevant gait events: toe-offs and heel-strikes). It enforces gait phase ordering in a cyclical manner: left swing, double-support 1, right swing, double support 2, left swing, ... and so on. 

The enforcing can be hard (no exceptions to cycle order) or soft (highly unlikely but not impossible exceptions). Soft enforcement can lead to slightly faster running times if accompanied by a sparse observation matrix (a sparse observation matrix is NOT recommended with hard enforcement: it may lead to ill-conditioned cases when subjects start/stop walking).

#to do:
Testing:
0) DEfine functions to convert from stance to events and viceversa, compute phase durations, and return summary stats of: phase duration and variability, 1-2-3-4 ordering exceptions (if any), force exceptions (e.g., non-zero force during swing, zero force during stance, if any).
1) Define metrics for comparison. Likelihood of data? Amount (and magnitude) of non-zero force samples during purported swing phases? Respecting the 1-2-3-4 phase ordering?
2) Test on non-thresholded (RAW) data.
3) Test on variety of subjects and trials.

Sys-id:
4) automate transition matrix estimation sys-id style


