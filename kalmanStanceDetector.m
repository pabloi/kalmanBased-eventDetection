function [pSmoothed,pUpdated,pPredicted,BW]=kalmanStanceDetector(force1,force2)
%This function implements a kalman-style filter (hidden markov chain state
%estimation) to infer stance from force traces.

BW=nanmean(force1+force2); %Median has better outlier rejection, but very skewed distributions may make it not a good estimator of weight (newton says mean(force in z direction) = weight because we are not floating away!
N=length(force1);
M=100; %Number of points where we estimate the probability distribution in the interval [-1 1]
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
sigma=.05;
pStateGivenPrevious=exp(-Dx.^2 ./ ((2*(1-abs(aux))+1)*sigma).^2); %Transition matrix
%pStateGivenPrevious=columnNormalize(pStateGivenPrevious);
pStateGivenPrevious=pStateGivenPrevious./sum(pStateGivenPrevious,1);
%Add higher chance of staying in place for -1 and 1 (point attractors)
w=.95; %This could be estimated from data. Should roughly be 2*sampling period/stride period.
pStateGivenPrevious(1,1)=w+(1-w)*pStateGivenPrevious(1,1);
pStateGivenPrevious(M,M)=w+(1-w)*pStateGivenPrevious(M,M);


%Forward pass (Kalman filter)
p0=ones(1,M)/M;
pPredicted=nan(N+1,M); %We can predict up to the Nth+1 sample
pPredicted(1,:)=p0;
pUpdated=nan(N,M);
for i=1:N
   %Update:
   pUpdated(i,:)=normalize(pObsGivenState(yObs(i,:),M).*pPredicted(i,:));
   %Predict:
   pPredicted(i+1,:)=sum(pStateGivenPrevious .* pUpdated(i,:),2)';
end

%Backward pass (Smoothing)
pSmoothed=nan(N,M);
pSmoothed(N,:)=pUpdated(N,:);
for i=(N-1):-1:1
    pSmoothed(i,:)=pUpdated(i,:) .* sum(pStateGivenPrevious .*(pSmoothed(i+1,:)./pPredicted(i+1,:))',1);
end

end

function p=normalize(p)
    p=p/sum(p(:));
end

function p=columnNormalize(p)
    p=p./sum(p,1);
end

function p=pObsGivenState(yObs,M)
%This is counter-intuitive: it computes the probability of an observation
%given a state, but does so as a function of the observation (which we
%have) and the unknown is the state
if any(yObs==0) %One of the measurements is 0: then the state is almost surely [1,0] or [0,1]
    aux=-sign(diff(yObs))-linspace(-1,1,M);
    p=1./(aux+1e-5);
else %None of the measurements is 0: use the difference to assign.
    aux=-diff(yObs)/abs(sum(yObs))-linspace(-1,1,M);
    sigma=.05;
    p=normalize(exp(-(aux.^2/sigma^2)) +1/(M*abs(diff(yObs))));
end
end