function [pSmoothed,pUpdated,pPredicted,pStateGivenPrevious,pObsGivenState,MAPsequence]=kalmanStanceDetectorv2(force1,force2)
%This function implements a kalman-style filter (hidden markov chain state
%estimation) to infer stance from force traces.

M=51; %Number of points where we estimate the probability distribution in the interval [-1 1]
D=101; %Number of points where we estimate the probability distribution of the output

%HMM chain description:
%The transition probability is: with probability=w we transition to a new
%state, with normal distribution around the current state, std=sigma and
%clipped to the acceptable interval [-1 1]. With prob=(1-w) we stay. w is 0
%for x=1/2 and 0.98 for x=+-1, with a quadratic shape.
aux=linspace(-1, 1,M);
Dx=aux' - aux;
sigma=.05;
pStateGivenPrevious=(exp(-Dx.^2 ./ ((2*(1-abs(aux))+1)*sigma).^2)); %Transition matrix
pStateGivenPrevious=(exp(-Dx.^2 ./ ((2*sigma).^2))); %Transition matrix
%Sparsify:
pStateGivenPrevious(pStateGivenPrevious<1e-5)=0;
pStateGivenPrevious=sparse(pStateGivenPrevious);
pStateGivenPrevious=columnNormalize(pStateGivenPrevious);

%Add higher chance of staying in place for -1 and 1 (point attractors)
w=.99; %This could be estimated from data. Should roughly be 2*sampling period/stride period.
pStateGivenPrevious(:,1)=[w;(1-w);zeros(M-2,1)];
pStateGivenPrevious(:,M)=[zeros(M-2,1);1-w;w];
%Also, make it hard to enter, as the probabilities are symmetric:
%pStateGivenPrevious(1,:)=pStateGivenPrevious(:,1)';
%pStateGivenPrevious(M,:)=pStateGivenPrevious(:,M)';

%Introduce histeresis:
L=tril(pStateGivenPrevious);
T=nan(2*M-2);
T(1:M,:)=[L(:,1:M-1),[L(end,1:M-1);zeros(M-1)]];
T(M:end,:)=[[L(end,1:M-1);zeros(M-2,M-1);],L(1:M-1,1:M-1)];
%T(sub2ind(size(T),1:2*M-3,2:2*M-2))=1e-20; %This makes the backward transition unlikely but not impossible. Comment for hard enforcement of cycles
T=sparse(T);
T=columnNormalize(T);

%Observation probability distribution:
measNoiseSigma=.03;
aux=[0:D-1]'/(D-1) - [0:M-1]/(M-1);
pObsGivenState=exp(-(aux.^2)/(2*measNoiseSigma^2))/sqrt(2*pi*measNoiseSigma^2);
%Sparsify (optional):
%pObsGivenState(pObsGivenState<1e-5)=0; %Reducing number of non-zero
%entries improves speed slightly. However, this is not recommended if there
%is hard-enforcement of the cycle, as it may lead to impossible state transitions.

pObsGivenState=columnNormalize(pObsGivenState);
pp=.05; %Prob of observing non-zero force in the single-stance phase
pObsGivenState(:,1)=[1-pp;pp/3;pp/3;pp/3;zeros(D-4,1)];
pObsGivenState(:,M)=[zeros(D-4,1);pp/3;pp/3;pp/3;1-pp];

O=[pObsGivenState(:,1:M-1),fliplr(pObsGivenState(:,2:M))];
%O=sparse(O); %If O is VERY sparse, this might be worth it in terms of
%speed, but probably not if O is more than 5% non-zero
O=columnNormalize(O);

p0=(ones(2*M-2,1)/(2*M-2));
observation=discretizeObs((force2-force1)./abs(force2+force1),D,[-1 1]);
[pPredicted, pUpdated, pSmoothed] = HMMstationaryInference(observation,O,T,p0);
[MAPsequence] = viterbi(observation,T,O,p0);
%Visualize matrices:
%figure; subplot(1,2,2); imagesc(T); title('Transition'); ylabel('Next state'); xlabel('Curr state'); subplot(1,2,1); imagesc(O); title('Observation');  ylabel('Obs'); xlabel('State');
end

function p=rowNormalize(p)
    %Normalization across columns
    p=p./sum(p,2);
end
function p=columnNormalize(p)
    %Normalization across columns
    p=p./sum(p,1);
end
