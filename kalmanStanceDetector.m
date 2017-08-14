function [pSmoothed,pUpdated,pPredicted,BW]=kalmanStanceDetector(force1,force2)
%This function implements a kalman-style filter (hidden markov chain state
%estimation) to infer stance from force traces.

BW=nanmean(force1+force2); %Median has better outlier rejection, but very skewed distributions may make it not a good estimator of weight (newton says mean(force in z direction) = weight because we are not floating away!
N=length(force1);
M=100; %Number of points where we estimate the probability distribution in the interval [-1 1]
D=100; %Number of points where we estimate the probability distribution of the output 
yObs=[force1 force2]/BW; %We have two observations that roughly will sum to 0, and both components are roughly always in [-1 1]
%yObs=(force1-force2)/BW;
%Alt idea: normalize the observations so that they sum to 1 on each point
%in time (i.e. we only care about force ratio)

%HMM chain description:
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

%Observation probability distribution:
measNoiseSigma=.03;
aux=[0:D-1]'/(D-1) - [0:M-1]/(M-1);
pObsGivenState=exp(-(aux.^2)/(2*measNoiseSigma^2))/sqrt(2*pi*measNoiseSigma^2);
pObsGivenState(:,1)= [.9; zeros(D-1,1)] + .1*pObsGivenState(:,1);%If state=1, then with high likelihood we observe 1
pObsGivenState(:,M)= [zeros(D-1,1);.9] + .1*pObsGivenState(:,M);%If state=1, then with high likelihood we observe 1





p0=ones(1,M)/M;
observation=-diff(yObs,[],2)./abs(sum(yObs,2));
%observation(yObs(1,:)==0)=1; %Assuming not both at 0 at the same time, ever
%observation(yObs(2,:)==0)=-1;
observation=round((D-1)*(observation+1)/2)+1; %Quantizing observations to [1,D] range
[pPredicted, pUpdated, pSmoothed] = genKFstationaryInference(observation,pObsGivenState,pStateGivenPrevious,p0);


figure; subplot(1,2,2); imagesc(pStateGivenPrevious); title('Transition'); ylabel('Next state'); xlabel('Curr state'); subplot(1,2,1); imagesc(pObsGivenState); title('Observation');  ylabel('Obs'); xlabel('State');
end

function p=normalize(p)
    p=p/sum(p(:));
end