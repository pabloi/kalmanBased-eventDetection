%% Add kalman engine to path
addpath(genpath('../matlab-sysID/')) %Change to wherever the matlab-linsys package is

%% Load
load ./data/C02forcesAndEventsTrial06.mat

%% Get force data from a single trial
Lfz=forces06(:,strcmp(Flabels,'LFz'));
Rfz=forces06(:,strcmp(Flabels,'RFz'));

%% Do the estimation
% [pS,pU,pP,T,O]=kalmanStanceDetector(Lfz,Rfz); %This actually estimates the % Fz force being exerted by the left leg (as opposed to the right)
% [~,MAPstate]=max(pS);
% MAPstate=MAPstate-1;
% swingL=MAPstate==0;
% swingR=MAPstate==100;
%% Alt estimation
[pS,pU,pP,T,O]=kalmanStanceDetectorv2(Lfz,Rfz); %This actually estimates the % Fz force being exerted by the left leg (as opposed to the right)
[~,MAPstate]=max(pS);
%[~,MAPstate]=max(pU);
MAPstate=2*(MAPstate-1);
MAPstate(MAPstate>100)=200-MAPstate(MAPstate>100);
swingL=MAPstate==0;
swingR=MAPstate==100;
%% Try EM:
%D=101;
%M=size(pS,2);
%observation=discretizeObs((Rfz-Lfz)./abs(Rfz+Lfz),D,[-1 1]);
%[O,T,states]=genEM(observation,ones(M,1)/M,O,T);
%[~,MAPstate]=max(states');
% MAPstate=MAPstate-1;
% swingL=MAPstate==0;
% swingR=MAPstate==100;
%% Compute associated events: LHS, RHS, LTO, RTO
eventsK=events06;
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
plot((MAPstate),'DisplayName','MAP estimate'); 
 linkaxes([p1 p2 p3 p4],'xy'); axis([1 size(Lfz,1) 0 100]); legend

%% Visualize aligned results
% sumfun=@mean;
% BW=median(Lfz+Rfz);
% figure; 
% p={pP',pU',pS'};
% tl={'Predicted','Updated','Smoothed'};
% ev=events06; %Use classic-detection events
% ev=eventsK; %Use kalma-filter based events -> If this is better, the forces should look better aligned to these events, both in mean/median and dispersion values
% for i=1:4 %Align to all 4 events
%     window=[-100:100];
%     aux=find(ev(:,i))+window;
%     for j=1:3
%     b=squeeze(mean(reshape(p{j}(aux(:),:),size(aux,1),size(aux,2),size(p{j},2)),1));
%     p1=subplot(4,4,i+(j-1)*4); imagesc(fliplr(b)'); ylabel('p(% on Right)');
%     set(p1,'XTick',[window(1),0:100:window(end)]-window(1)+1,'XTickLabel',{num2str(window(1)),Elabels{i}(end-2:end),'+100','+200'})
%     if j==1
%         title(['Aligned to ' Elabels{i}(end-2:end)])
%     end
%     end
%     p1=subplot(4,4,i+12);     %j=4, plot the actual data:
%     hold on; 
%     l1=plot(100*sumfun(Lfz(aux),1)/BW,'DisplayName','% BW on left','LineWidth',2); 
%     l2=plot(100-(100*sumfun(Rfz(aux),1)/BW),'DisplayName','100-%BW on right','LineWidth',2);
%     l3=plot(((sumfun(MAPstate(aux),1))),'DisplayName','MAP estimate','LineWidth',2); 
%     if mod(i,2)==1
%     plot(100*(Lfz(aux))'/BW,'DisplayName','% BW on left','Color',l1.Color); 
%     axis([0 200 0 30])
%     else
%     plot(100-(100*(Rfz(aux))'/BW),'DisplayName','100-%BW on right','Color',l2.Color);
%     axis([0 200 70 100])
%     end
%     plot((((MAPstate(aux)))'),'DisplayName','MAP estimate','Color',l3.Color); 
%     %TODO add some sense of dispersion to the plots
%     set(p1,'XTick',[window(1),0:100:window(end)]-window(1)+1,'XTickLabel',{num2str(window(1)),Elabels{i}(end-2:end),'+100','+200'})
%     grid on
%     if i==4
%         legend([l1,l2,l3],'Location','SouthEast')
%     end
% end

%% Visualize aligned results to the two event types
sumfun=@median;
dispfun=@iqr;
dispfunD=@(x) median(x)-prctile(x,16,1);
dispfunU=@(x) prctile(x,84,1)-median(x);
dispfunD2=@(x) median(x)-prctile(x,97.5,1);
dispfunU2=@(x) prctile(x,2.5,1)-median(x);
sumfun=@mean;
%dispfunD=@std;
%dispfunU=@std;
BW=median(Lfz+Rfz);
figure; 
for j=1:2
    switch j
        case 1
            ev=events06;
            en='Classic';
        case 2
            ev=eventsK;
            en='Kalman';
    end
    window=[-900:900];
    for i=1:4 %Align to all 4 events
        LHS=find(ev([(-window(1)+1):end-window(end)],i))+window-window(1)-1; %MAking sure we have enough data to plot
        for k=1:3 %Lfz, Rfz, state
            switch k
                case 1
                    data=100*Lfz(LHS)/BW; %Scaling for viz
                    nn='% BW on left';
                case 2
                    data=100-100*Rfz(LHS)/BW;
                    nn='100-%BW on right';
                case 3
                    data=MAPstate(LHS);
                    nn='MAP estimate';
            end
            p1=subplot(3,4,i+(k-1)*4);     %j=4, plot the actual data:
            hold on; 
           l=plot(sumfun(data),'DisplayName',en,'LineWidth',2);
           patch([1:size(data,2),size(data,2):-1:1],[sumfun(data)+dispfunU2(data),fliplr(sumfun(data)-dispfunD2(data))],l.Color,'FaceAlpha',.1,'EdgeColor','none','DisplayName','90% CI (\pm 2 \sigma)')
           patch([1:size(data,2),size(data,2):-1:1],[sumfun(data)+dispfunU(data),fliplr(sumfun(data)-dispfunD(data))],l.Color,'FaceAlpha',.3,'EdgeColor','none','DisplayName','68% CI (\pm 1 \sigma)')
           if j==1
              %plot(data','Color',l.Color) 
           end
           if i==1
               ylabel(nn)
           end
                   set(p1,'XTick',[window(1),-300:300:window(end)]-window(1)+1,'XTickLabel',{num2str(window(1)),'-300',Elabels{i}(end-2:end),'+300','+600'})
        grid on
        axis tight
        end

        if i==4
            legend('Location','SouthEast')
        end
    end

end
%% Visualize distribution of all phases of gait:
figure
for j=1:2
    switch j
        case 1
            ev=events06;
            en='Classic';
        case 2
            ev=eventsK;
            en='Kalman';
    end
    LHS=find(ev(:,1)); %Starting and ending LHS
    RHS=nan(size(LHS));
    RTO=nan(size(LHS));
    LTO=nan(size(LHS));
    DS1=nan(size(LHS));
    DS2=nan(size(LHS));
    sstL=nan(size(LHS));
    sstR=nan(size(LHS));
   for k=1:numel(LHS)-1
       RHS(k)=find(ev(LHS(k):LHS(k+1),2),1,'first')+LHS(k);
       RTO(k)=find(ev(LHS(k):LHS(k+1),4),1,'first')+LHS(k);
       LTO(k)=find(ev(LHS(k):LHS(k+1),3),1,'first')+LHS(k);
       DS1(k)=RTO(k)-LHS(k);
       DS2(k)=LTO(k)-RHS(k);
       sstR(k)=LHS(k+1)-LTO(k)-1;
       sstL(k)=RHS(k)-RTO(k);
   end
   for k=1:4 %4 phases of gait, draw hist
       switch k
           case 1
               data=DS1;
               nn='LHS to RTO';
           case 2
               data=DS2;
               nn='RHS to LTO';
           case 3
               data=sstR;
               nn='Single Stance R';
           case 4
               data=sstL;
               nn='Single Stance L';
       end
      subplot(2,2,k)
      cc=get(gca,'ColorOrder');
      hold on
      h=histogram(data);
      title(nn)
      text(1.1*nanmean(data),30-2*j,['\mu = ' num2str(nanmean(data)) ' \pm ' num2str(nanstd(data))],'Color',cc(j,:))
   end
end

%%  Print difference in forces during 'swing'
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

