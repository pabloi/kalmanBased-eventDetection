%% Add kalman engine to path
addpath(genpath('../matlab-sysID/')) %Change to wherever the matlab-linsys package is

%% Load
load ./data/C02forcesAndEventsTrial06.mat

%% Get force data from a single trial
Lfz=forces06(:,strcmp(Flabels,'LFz'));
Rfz=forces06(:,strcmp(Flabels,'RFz'));

%% Do the estimation
[pS,pU,pP,BW]=kalmanStanceDetector(Lfz,Rfz); %This actually estimates the % Fz force being exerted by the left leg (as opposed to the right)
[~,MAPstate]=max(pS');
%% Compute associated events: LHS, RHS, LTO, RTO
eventsK=events06;
eventsK(2:end,1)= MAPstate(2:end)>1 & MAPstate(1:end-1)==1; %LHS
eventsK(2:end,2)= MAPstate(2:end)<100 & MAPstate(1:end-1)==100; %RHS
eventsK(2:end,3)= MAPstate(2:end)==1 & MAPstate(1:end-1)>1; %LTO
eventsK(2:end,4)= MAPstate(2:end)==100 & MAPstate(1:end-1)<100; %LTO
%% Visualize results
figure; 
p1=subplot(4,1,1); imagesc(pP'); title('Predicted State')
p2=subplot(4,1,2); imagesc(pU'); title('Updated State')
p3=subplot(4,1,3); imagesc(pS'); title('Smoothed State')
p4=subplot(4,1,4); hold on; plot((100*Lfz/BW),'DisplayName','% BW on left'); 
plot(100-(100*Rfz/BW),'DisplayName','100-%BW on right');
plot(((MAPstate-1)*100/99),'DisplayName','MAP estimate'); 
linkaxes([p1 p2 p3 p4],'xy'); axis([1 size(Lfz,1) 0 100]); legend

%% Visualize aligned results
figure; 
p={pP,pU,pS};
tl={'Predicted','Updated','Smoothed'};
ev=events06; %Use classic-detection events
ev=eventsK; %Use kalma-filter based events -> If this is better, the forces should look better aligned to these events, both in mean/median and dispersion values
for i=1:4 %Align to all 4 events
    window=[-30:30];
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
    plot(100*(Lfz(aux))'/BW,'DisplayName','% BW on left','Color',l1.Color); 
    plot(100-(100*(Rfz(aux))'/BW),'DisplayName','100-%BW on right','Color',l2.Color);
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

%% To DO: visualize actual forces to distribution of estimated forces.-> Is the actual force close to the MAP?