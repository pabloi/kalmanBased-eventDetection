function [pSmoothed,pUpdated,pPredicted,BW]=percentFZsmootherv2(force1,force2)
%Kalman-style smoothing of data to infer % contribution of each of two
%forces, AND ITS DERIVATIVE, to the net force. 
%Useful to infer denoised stance from force traces. 
%The kalman engine works on discrete 1D states and 1D observations only, 
%so the problem needs to be set up as such.
%See also: percentFZsmoother, genKFstationaryInference
error('Unimplemented')

end