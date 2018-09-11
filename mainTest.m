%% Add kalman engine to path
addpath(genpath('../matlab-sysID/')) %Change to wherever the matlab-linsys package is

%% Load
load ./data/C02forcesAndEventsTrial06.mat

%% Get force data from a single trial
Lfz=forces06(:,strcmp(Flabels,'LFz'));
Rfz=forces06(:,strcmp(Flabels,'RFz'));

%% Do the estimation
[pS,pU,pP,T,O]=kalmanStanceDetector(Lfz,Rfz); %This actually estimates the % Fz force being exerted by the left leg (as opposed to the right)
[~,MAPstate]=max(pS');
%% Try EM:
D=100;
M=size(pS,2);
observation=discretizeObs((Rfz-Lfz)./abs(Rfz+Lfz),D,[-1 1]);
[O,T]=genEM(observation,ones(M,1)/M,O,T);
%% Compute associated events: LHS, RHS, LTO, RTO
eventsK=events06;
eventsK(2:end,1)= MAPstate(2:end)>1 & MAPstate(1:end-1)==1; %LHS
eventsK(2:end,2)= MAPstate(2:end)<100 & MAPstate(1:end-1)==100; %RHS
eventsK(2:end,3)= MAPstate(2:end)==1 & MAPstate(1:end-1)>1; %LTO
eventsK(2:end,4)= MAPstate(2:end)==100 & MAPstate(1:end-1)<100; %RTO
%% Visualize results
% figure; 
% p1=subplot(4,1,1); imagesc(pP'); title('Predicted State')
% p2=subplot(4,1,2); imagesc(pU'); title('Updated State')
% p3=subplot(4,1,3); imagesc(pS'); title('Smoothed State')
% p4=subplot(4,1,4); hold on; plot((100*Lfz/BW),'DisplayName','% BW on left'); 
% plot(100-(100*Rfz/BW),'DisplayName','100-%BW on right');
% plot(((MAPstate-1)*100/99),'DisplayName','MAP estimate'); 
% linkaxes([p1 p2 p3 p4],'xy'); axis([1 size(Lfz,1) 0 100]); legend

%% Visualize aligned results
figure; 
p={pP,pU,pS};
tl={'Predicted','Updated','Smoothed'};
ev=events06; %Use classic-detection events
ev=eventsK; %Use kalma-filter based events -> If this is better, the forces should look better aligned to these events, both in mean/median and dispersion values
for i=1:4 %Align to all 4 events
    window=[-200:200];
    aux=find(ev(:,i))+window;
    for j=1:3
    b=squeeze(mean(reshape(p{j}(aux(:),:),size(aux,1),size(aux,2),size(p{j},2)),1));
    p1=subplot(4,4,i+(j-1)*4); imagesc(fliplr(b)'); ylabel('p(% on Right)');
    set(p1,'XTick',[window(1),0:100:window(end)]-window(1)+1,'XTickLabel',{num2str(window(1)),Elabels{i}(end-2:end),'+100','+200'})
    if j==1
        title(['Aligned to ' Elabels{i}(end-2:end)])
    end
    end
    %j=4, plot the data:
    p1=subplot(4,4,i+12);
    hold on; 
    l1=plot(100*mean(Lfz(aux),1)/BW,'DisplayName','% BW on left','LineWidth',2); 
    l2=plot(100-(100*mean(Rfz(aux),1)/BW),'DisplayName','100-%BW on right','LineWidth',2);
    l3=plot(((mean(MAPstate(aux),1)-1)*100/99),'DisplayName','MAP estimate','LineWidth',2); 
    %plot(100*median(Lfz(aux),1)/BW,'DisplayName','median % BW on left'); 
    %plot(100-(100*median(Rfz(aux),1)/BW),'DisplayName','median 100-%BW on right');
    %plot(((median(MAPstate(aux),1)-1)*100/99),'DisplayName','median MAP estimate');
    if mod(i,2)==1
    plot(100*(Lfz(aux))'/BW,'DisplayName','% BW on left','Color',l1.Color); 
    axis([100 300 0 15])
    else
    plot(100-(100*(Rfz(aux))'/BW),'DisplayName','100-%BW on right','Color',l2.Color);
    axis([100 300 85 100])
    end
    plot((((MAPstate(aux))-1)'*100/99),'DisplayName','MAP estimate','Color',l3.Color); 
    %TODO add some sense of dispersion to the plots
    set(p1,'XTick',[window(1),0:100:window(end)]-window(1)+1,'XTickLabel',{num2str(window(1)),Elabels{i}(end-2:end),'+100','+200'})
    grid on
    if i==4
        legend([l1,l2,l3],'Location','SouthEast')
    end
end

%p4=subplot(4,1,4); hold on; plot((100*Lfz/BW),'DisplayName','% BW on left'); 
%plot(100-(100*Rfz/BW),'DisplayName','100-%BW on right');
%plot(((MAPstate-1)*100/99),'DisplayName','MAP estimate'); 
%linkaxes([p1 p2 p3 p4],'xy'); axis([1 size(Lfz,1) 0 100]); 

%%  Visualize difference in forces during 'stance'
ev=full(events06);
stanceL{1}=1+cumsum(ev(:,1))-cumsum(ev(:,3));
stanceR{1}=1+cumsum(ev(:,2))-cumsum(ev(:,4));
DS_LR{1}=1+cumsum(ev(:,1))-cumsum(ev(:,4));
DS_RL{1}=1+cumsum(ev(:,2))-cumsum(ev(:,3));
ev=full(eventsK);
stanceL{2}=cumsum(ev(:,1))-cumsum(ev(:,3));
stanceL{2}=stanceL{2}-median(stanceL{2})+1;
stanceR{2}=cumsum(ev(:,2))-cumsum(ev(:,4));
stanceR{2}=stanceR{2}-median(stanceR{2})+1;
DS_LR{2}=1+cumsum(ev(:,1))-cumsum(ev(:,4));
DS_RL{2}=1+cumsum(ev(:,2))-cumsum(ev(:,3));

%Non-zero forces measured during alleged 'swing': (Type I error), sometimes happens
sum(Lfz(stanceL{1}==0)~=0)
sum(Lfz(stanceL{2}==0)~=0)
sum(Rfz(stanceR{1}==0)~=0)
sum(Rfz(stanceR{2}==0)~=0)

%0-force samples during alleged 'stance': (Type II error), this should almost never happen
sum(Lfz(stanceL{1}==1)==0)
sum(Lfz(stanceL{2}==1)==0)
sum(Rfz(stanceR{1}==1)==0)
sum(Rfz(stanceR{2}==1)==0)

