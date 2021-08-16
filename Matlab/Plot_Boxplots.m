clear all

figure(1); clf; CL = get(gca,'colororder'); FS = 15; set(gcf,'color','w');

load TEMP Impact_Large MV
Impact = Impact_Large;
NumSpp = size(Impact,2);
SP = subplot('position',[0.1 0.58 0.5 0.35]); hold on; d = 0.1;
Q = quantile(Impact,[0.16 0.5 0.84],3);
REF = mean(mean(mean(Q(:,1,:))));
Q = Q./REF;
for m = 1:length(MV)
    for s = 1:NumSpp
        ppp(m) = plot(m+5*s+[-d d nan 0 0 nan -d d],squeeze(Q(m,s,[1 1 1 1 3 3 3 3])),'k','linewidth',2,'color',CL(7-m,:));
        plot(m+5*s,squeeze(Q(m,s,2)),'k.','markersize',12,'color',CL(7-m,:))
    end
end

set(gca,'xtick',[7 12 17 22],'xticklabel',{'Global','Long MDD','Short MDD','Brooder'},'fontsize',FS);
L = legend(ppp,'High correlation disturbances','Low correlation disturbances','No correlation');
set(L,'location','northwest','box','off','fontsize',FS)
title('Large magnitude disturbance','fontsize',FS+2,'fontweight','normal')
ylabel('Relative impact','fontsize',FS+2,'fontweight','normal')
ylim([0.9 1.8])

clearvars -except CL FS
load TEMP Impact_Small MV
Impact = Impact_Small;
NumSpp = size(Impact,2);

SP = subplot('position',[0.1 0.07 0.5 0.35]); hold on; d = 0.1;
Q = quantile(Impact,[0.16 0.5 0.84],3);
REF = mean(mean(mean(Q(:,1,:))));
Q = Q./REF;
for m = 1:length(MV)
    for s = 1:NumSpp
        ppp(m) = plot(m+5*s+[-d d nan 0 0 nan -d d],squeeze(Q(m,s,[1 1 1 1 3 3 3 3])),'k','linewidth',2,'color',CL(7-m,:));
        plot(m+5*s,squeeze(Q(m,s,2)),'k.','markersize',12,'color',CL(7-m,:))
    end
end

set(gca,'xtick',[7 12 17 22],'xticklabel',{'Global','Long MDD','Short MDD','Brooder'},'fontsize',FS);
L = legend(ppp,'High correlation disturbances','Low correlation disturbances','No correlation');
set(L,'location','northwest','box','off','fontsize',FS)
title('Small magnitude disturbance','fontsize',FS+2,'fontweight','normal')
ylabel('Relative impact','fontsize',FS+2,'fontweight','normal')
ylim([0.9 2.8])

SP = subplot('position',[0.63 0 0.35 0.99]); hold on; MS = 6; axis off
load crcb_domain_full centroid

Small = 0;
if Small == 1
    F = find(centroid(:,2) < -16);
    centroid(F,:) = [];
end
NumReefs = length(centroid);
Distance = squareform(pdist(centroid))*110;

% Pick a central reef
CR = 700;
% CR = randsample(NumReefs,1)
[~,I] = sort(Distance(CR,:),'ascend');

SHF = 16;
XSHF = 5;

% Small, high
M = [0 -1];
Spacing = 1; 
NumDamaged = 250;
Damaged = I(1:Spacing:Spacing*NumDamaged);
plot(XSHF*M(1)+centroid(:,1),      SHF*M(2)+centroid(:,2),'.','markersize',MS,'color',0.7.*ones(1,3))
plot(XSHF*M(1)+centroid(Damaged,1),SHF*M(2)+centroid(Damaged,2),'.','markersize',MS+1,'color',CL(6,:))

% Small, low
M = [1 -1];
Spacing = 3; 
Damaged = I(1:Spacing:Spacing*NumDamaged);
plot(XSHF*M(1)+centroid(:,1),      SHF*M(2)+centroid(:,2),'.','markersize',MS,'color',0.7.*ones(1,3))
plot(XSHF*M(1)+centroid(Damaged,1),SHF*M(2)+centroid(Damaged,2),'.','markersize',MS+1,'color',CL(5,:))

% Small, none
M = [2 -1];
Damaged = randsample(NumReefs,NumDamaged,0);
plot(XSHF*M(1)+centroid(:,1),      SHF*M(2)+centroid(:,2),'.','markersize',MS,'color',0.7.*ones(1,3))
plot(XSHF*M(1)+centroid(Damaged,1),SHF*M(2)+centroid(Damaged,2),'.','markersize',MS+1,'color',CL(4,:))

% large, high
M = [0 0];
Spacing = 1; 
NumDamaged = 500;
Damaged = I(1:Spacing:Spacing*NumDamaged);
plot(XSHF*M(1)+centroid(:,1),      SHF*M(2)+centroid(:,2),'.','markersize',MS,'color',0.7.*ones(1,3))
plot(XSHF*M(1)+centroid(Damaged,1),SHF*M(2)+centroid(Damaged,2),'.','markersize',MS+1,'color',CL(6,:))

% large, low
M = [1 0];
Spacing = 4; 
Damaged = I(1:Spacing:Spacing*NumDamaged);
plot(XSHF*M(1)+centroid(:,1),      SHF*M(2)+centroid(:,2),'.','markersize',MS,'color',0.7.*ones(1,3))
plot(XSHF*M(1)+centroid(Damaged,1),SHF*M(2)+centroid(Damaged,2),'.','markersize',MS+1,'color',CL(5,:))

% Large, none
M = [2 0];
Damaged = randsample(NumReefs,NumDamaged,0);
plot(XSHF*M(1)+centroid(:,1),      SHF*M(2)+centroid(:,2),'.','markersize',MS,'color',0.7.*ones(1,3))
plot(XSHF*M(1)+centroid(Damaged,1),SHF*M(2)+centroid(Damaged,2),'.','markersize',MS+1,'color',CL(4,:))

ylim([-42 -8]+0.2)
% xlim([141 151])














