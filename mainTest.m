%% Load
load ../../PhD/lab/rawData/synergies/mat/C0001.mat

%% Get force data from a single trial
trial=5;
ff=expData.data{trial}.GRFData.getDataAsVector({'LFz','RFz'});

%% Do the estimation
[pS,pU,pP,BW]=kalmanStanceDetector(ff(:,1),ff(:,2));
[~,MAPstate]=max(pS');

%% Visualize results
figure; p1=subplot(4,1,1); imagesc(pP'); p2=subplot(4,1,2); imagesc(pU'); p3=subplot(4,1,3); imagesc(pS'); p4=subplot(4,1,4); hold on; plot((100*ff(:,1)/BW)); plot(100-(100*ff(:,2)/BW)); plot(((MAPstate-1)*100/99)); linkaxes([p1 p2 p3 p4],'xy'); axis([1 size(ff,1) 0 100])

%% To DO: visualize actual forces to distribution of estimated forces.-> Is the actual force close to the MAP?