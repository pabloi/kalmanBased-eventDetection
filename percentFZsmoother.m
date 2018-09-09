function [pSmoothed,pUpdated,pPredicted,BW]=percentFZsmoother(force1,force2)
%Kalman-style smoothing of data to infer % contribution of each of two
%forces to the net force. Useful to infer denoised stance from force traces. 
%The kalman engine works on discrete 1D states and 1D observations only, 
%so the problem needs to be set up as such.
%See also: percentFZsmootherv2, genKFstationaryInference

BW=nanmean(force1+force2); %Median has better outlier rejection, but very skewed distributions may make it not a good estimator of weight (newton says mean(force in z direction) = weight because we are not floating away!
N=length(force1);
M=100; %Number of points where we estimate the probability distribution in the interval [-1 1]
D=100; %Number of points where we estimate the probability distribution of the output (% Fz)

%HMM chain description: DEFINE STATE TRANSITION DENSITY
%The transition probability is: with probability=w we transition to a new
%state, with normal distribution around the current state, std=sigma and
%clipped to the acceptable interval [-1 1]. With prob=(1-w) we stay. w is 0
%for x=1/2 and 0.98 for x=+-1, with a quadratic shape.
aux=linspace(-1, 1,M);
Dx=aux' - aux;
sigma=.01;
pStateGivenPrevious=exp(-Dx.^2 ./ ((2*(1-abs(aux))+1)*sigma).^2); %Transition matrix
%pStateGivenPrevious=columnNormalize(pStateGivenPrevious);
pStateGivenPrevious=pStateGivenPrevious./sum(pStateGivenPrevious,1);
%Add higher chance of staying in place for -1 and 1 (point attractors)
w=.9; %This could be estimated from data. Should roughly be 2*sampling period/stride period.
pStateGivenPrevious(1,1)=w+(1-w)*pStateGivenPrevious(1,1);
pStateGivenPrevious(M,M)=w+(1-w)*pStateGivenPrevious(M,M);

%Observation probability distribution: DEFINE OBSERVATION DENSITY
measNoiseSigma=.03;
aux=[0:D-1]'/(D-1) - [0:M-1]/(M-1);
pObsGivenState=exp(-(aux.^2)/(2*measNoiseSigma^2))/sqrt(2*pi*measNoiseSigma^2);
pObsGivenState(:,1)= [.9; zeros(D-1,1)] + .1*pObsGivenState(:,1);%If state=1, then with high likelihood we observe 1
pObsGivenState(:,M)= [zeros(D-1,1);.9] + .1*pObsGivenState(:,M);%If state=1, then with high likelihood we observe 1

%QUANTIZE OBSERVATIONS
yObs=[force1 force2]/BW;
p0=ones(1,M)/M;
observation=-diff(yObs,[],2)./abs(sum(yObs,2)); % %Fz observed on first force, always bound to [-1 1]
observation=round((D-1)*(observation+1)/2)+1; %Quantizing observations to [1,D] range
[pPredicted, pUpdated, pSmoothed] = genKFstationaryInference(observation,pObsGivenState,pStateGivenPrevious,p0);

%Visualize transition and observation matrices, if desired:
%figure; subplot(1,2,2); imagesc(pStateGivenPrevious); title('Transition'); ylabel('Next state'); xlabel('Curr state'); subplot(1,2,1); imagesc(pObsGivenState); title('Observation');  ylabel('Obs'); xlabel('State');
end

function p=normalize(p)
    p=p/sum(p(:));
end