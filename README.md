kalmanBased-eventDetection

#to do:
Testing:
Define metrics for comparison. Likelihood of data? Amount (and magnitude) of non-zero force samples during purported swing phases? Respecting the 1-2-3-4 phase ordering?

Given transition matrices:
1) stance/swing detection from ipsilateral z-force readings
2) singleL/DS1/singleR/DS2 detection from bilateral z-force readings
3) singleL/DS1/singleR/DS2 detection from bilateral 3D force readings
4) 2-stage filter: first filter z-force readings to estimate current value and derivative, then use the filtered estimate to repeat 2)

Sys-id:
4) automate transition matrix estimation sys-id style
