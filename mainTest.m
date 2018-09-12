%% Add kalman engine to path
addpath(genpath('../matlab-sysID/')) %Change to wherever the matlab-linsys package is

%% Load
load ./data/C02forcesAndEventsTrial06.mat

%% Get force data from a single trial
Lfz=forces06(:,strcmp(Flabels,'LFz'));
Rfz=forces06(:,strcmp(Flabels,'RFz'));

%% Do the estimation
[pS,pU,pP,T,O]=kalmanStanceDetector(Lfz,Rfz); %This actually estimates the % Fz force being exerted by the left leg (as opposed to the right)
[~,MAPstate]=max(pS);
MAPstate=MAPstate-1;
swingL=MAPstate==0;
swingR=MAPstate==100;
%% Alt estimation
% [pS,pU,pP,T,O]=kalmanStanceDetectorv2(Lfz,Rfz); %This actually estimates the % Fz force being exerted by the left leg (as opposed to the right)
% [~,MAPstate]=max(pS);
% %[~,MAPstate]=max(pU);
% MAPstate=2*(MAPstate-1);
% MAPstate(MAPstate>100)=200-MAPstate(MAPstate>100);
% swingL=MAPstate==0;
% swingR=MAPstate==100;
%% Try EM:
%D=100;
%M=size(pS,2);
%observation=discretizeObs((Rfz-Lfz)./abs(Rfz+Lfz),D,[-1 1]);
%[O,T,states]=genEM(observation,ones(M,1)/M,O,T);
%[~,MAPstate]=max(states');
%% Compute associated events: LHS, RHS, LTO, RTO
eventsK=events06;
% eventsK(2:end,1)= MAPstate(2:end)>1 & MAPstate(1:end-1)==1; %LHS
% eventsK(2:end,2)= MAPstate(2:end)<100 & MAPstate(1:end-1)==100; %RHS
% eventsK(2:end,3)= MAPstate(2:end)==1 & MAPstate(1:end-1)>1; %LTO
% eventsK(2:end,4)= MAPstate(2:end)==100 & MAPstate(1:end-1)<100; %RTO
eventsK(2:end,1)= ~swingL(2:end) & swingL(1:end-1); %LHS
eventsK(2:end,2)= ~swingR(2:end) & swingR(1:end-1); %RHS
eventsK(2:end,3)= swingL(2:end) & ~swingL(1:end-1); %LTO
eventsK(2:end,4)= swingR(2:end) & ~swingR(1:end-1); %RTO
%% Visualize results
figure; 
p1=subplot(4,1,1); imagesc(pP); title('Predicted State')
p2=subplot(4,1,2); imagesc(pU); title('Updated State')
p3=subplot(4,1,3); imagesc(pS); title('Smoothed State')
p4=subplot(4,1,4); hold on; plot((100*Lfz/BW),'DisplayName','% BW on left'); 
plot(100-(100*Rfz/BW),'DisplayName','100-%BW on right');
plot(((MAPstate-1)*100/99),'DisplayName','MAP estimate'); 
 linkaxes([p1 p2 p3 p4],'xy'); axis([1 size(Lfz,1) 0 100]); legend

%% Visualize aligned results
BW=median(Lfz+Rfz);
figure; 
p={pP',pU',pS'};
tl={'Predicted','Updated','Smoothed'};
ev=events06; %Use classic-detection events
ev=eventsK; %Use kalma-filter based events -> If this is better, the forces should look better aligned to these events, both in mean/median and dispersion values
for i=1:4 %Align to all 4 events
    window=[-100:100];
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
    l3=plot(((mean(MAPstate(aux),1))),'DisplayName','MAP estimate','LineWidth',2); 
    %plot(100*median(Lfz(aux),1)/BW,'DisplayName','median % BW on left'); 
    %plot(100-(100*median(Rfz(aux),1)/BW),'DisplayName','median 100-%BW on right');
    %plot(((median(MAPstate(aux),1)-1)*100/99),'DisplayName','median MAP estimate');
    if mod(i,2)==1
    plot(100*(Lfz(aux))'/BW,'DisplayName','% BW on left','Color',l1.Color); 
    axis([0 200 0 30])
    else
    plot(100-(100*(Rfz(aux))'/BW),'DisplayName','100-%BW on right','Color',l2.Color);
    axis([0 200 70 100])
    end
    plot((((MAPstate(aux)))'),'DisplayName','MAP estimate','Color',l3.Color); 
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

%%  Visualize difference in forces during 'swing'
ev=full(events06);
swL{1}=cumsum(ev(:,3))-cumsum(ev(:,1));
swR{1}=cumsum(ev(:,4))-cumsum(ev(:,2));
ev=full(eventsK);
swL{2}=swingL;
swR{2}=swingR;

%Non-zero forces measured during alleged 'swing': (Type I error), sometimes happens
sum(Lfz(swL{1}==1)~=0)
sum(Lfz(swL{2}==1)~=0)
sum(Rfz(swR{1}==1)~=0)
sum(Rfz(swR{2}==1)~=0)

%0-force samples during alleged 'stance': (Type II error), this should almost never happen
sum(Lfz(swL{1}==0)==0)
sum(Lfz(swL{2}==0)==0)
sum(Rfz(swR{1}==0)==0)
sum(Rfz(swR{2}==0)==0)

