function [pSmoothed,pUpdated,pPredicted,pStateGivenPrevious,pObsGivenState]=kalmanStanceDetector(force1,force2)
%This function implements a kalman-style filter (hidden markov chain state
%estimation) to infer stance from force traces.

M=100; %Number of points where we estimate the probability distribution in the interval [-1 1]
D=100; %Number of points where we estimate the probability distribution of the output 

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
pp=.97;
pObsGivenState(:,1)= [pp; zeros(D-1,1)] + (1-pp)*pObsGivenState(:,1);%If state=1, then with high likelihood we observe 1
pObsGivenState(:,M)= [zeros(D-1,1);pp] + (1-pp)*pObsGivenState(:,M);%If state=1, then with high likelihood we observe 1





p0=ones(1,M)/M;
observation=discretizeObs((force2-force1)./abs(force2+force1),D,[-1 1]);
[pPredicted, pUpdated, pSmoothed] = genKFstationaryInference(observation,pObsGivenState,pStateGivenPrevious,p0);
%Visualize matrices:
%figure; subplot(1,2,2); imagesc(pStateGivenPrevious); title('Transition'); ylabel('Next state'); xlabel('Curr state'); subplot(1,2,1); imagesc(pObsGivenState); title('Observation');  ylabel('Obs'); xlabel('State');
end