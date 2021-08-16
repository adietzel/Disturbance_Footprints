clear all

figure(1); clf; CL = get(gca,'colororder'); FS = 15; set(gcf,'color','w');

% Load the locations of the GBR centroids
load crcb_domain_full centroid
return
Small = 0;
if Small == 1
    F = find(centroid(:,2) < -18);
    centroid(F,:) = [];
end

F = find(isnan(centroid(:,1)));
centroid(F,:) = [];
Distance = squareform(pdist(centroid))*110;
NumReefs = length(centroid)
NumSpp = 4;
NumReps = 3;
f = 10e-2;
MV = [1 3 nan]; % Spacing between reefs to mimic spatial autocorrelation

% Load the best fit kernels
load FittingResults p_best L
% SP = [1478 1472]-04;
SP = [792 795];
k = [-1e-6; p_best(SP,2); -0.9];

for i = 1:NumSpp
    C(:,:,i) = exp(k(i).*Distance);
    
    %     % Normalise by the entire matrix
    %     C(:,:,i) = C(:,:,i)./sum(sum(C(:,:,i)));
    % Normalise by each row
    C(:,:,i) = C(:,:,i)./repmat(sum(C(:,:,i),2),1,NumReefs);
end

%% === LARGE MAGNITUDE DISTURBANCE ===
NumDamaged = 500;
for m = 1:length(MV)
    MaxSpace = MV(m);
    for reps = 1:NumReps
        disp(['Large disturbance, spacing = ' num2str(MV(m)) '; rep # ' num2str(reps)])
        
        % Pick a central reef
        CR = randsample(NumReefs,1);
        [~,I] = sort(Distance(CR,:),'ascend');
        
        % Damage some reefs
        Spacing(1,reps) = MaxSpace;%randi(MaxSpace);
        if MaxSpace == nan
            Damaged = randsample(NumReefs,NumDamaged,0);
        else
            Damaged = I(2:Spacing(1,reps):Spacing(1,reps)*NumDamaged);
        end
        
        % Let them recover
        N = ones(NumSpp,NumReefs);
        for t = 1:1e4
            
            if t == 5
                N(:,Damaged) = 0.4;
            end
            
            for s = 1:NumSpp
                
                % Local-growth driven recovery at a 5% rate
                N(s,:) = min(1,N(s,:)*1.05);
                
                % Background mortality
                Mortrate = ones(1,NumReefs);
                Mortrate(rand(size(Mortrate))<0.4) = 0.85;
                N(s,:) = N(s,:).*Mortrate;
                
                % Quadratic fertilisation success for local production
                LarvalOutput = f.*N(s,:).^2;
                
                % Dispersal-driven recovery
                N(s,:) = min(1,N(s,:) + (LarvalOutput*C(:,:,s)).*(1-N(s,:)));
                SN(s,t) = sum(N(s,:));
            end
            if t > 5 & sum(abs(SN(:,end)-SN(:,1))) < 1
                break
            end
        end
        for s = 1:NumSpp
            Impact(m,s,reps) = sum(abs(SN(s,1) - SN(s,:)));
        end
        
        clearvars -except reps centroid C Num* Distance CL Impact Spacing MaxSpace k m MV f FS
    end
end
subplot(2,1,1); hold on; d = 0.1;
Q = quantile(Impact,[0.16 0.5 0.84],3)
for m = 1:length(MV)
    for s = 1:NumSpp
        ppp(m) = plot(m+5*s+[-d d nan 0 0 nan -d d],squeeze(Q(m,s,[1 1 1 1 3 3 3 3])),'k','linewidth',2,'color',CL(7-m,:));
        plot(m+5*s,squeeze(Q(m,s,2)),'k.','markersize',12,'color',CL(7-m,:))
    end
end

set(gca,'xtick',[7 12 17 22],'xticklabel',{'Global disperser','Long PLD spawner','Short PLD spawner','Brooder'},'fontsize',FS);
L = legend(ppp,'High correlation disturbances','Low correlation disturbances','No correlation');
set(L,'location','northwest','box','off','fontsize',FS)
title('Large magnitude disturbance','fontsize',FS+2,'fontweight','normal')

Impact_Large = Impact;
save TEMP Impact* MV centroid

%% === SMALL MAGNITUDE DISTURBANCE ===
clearvars -except centroid C Num* Distance CL k MV f FS
NumDamaged = 100;

for m = 1:length(MV)
    MaxSpace = MV(m);
    for reps = 1:NumReps
        disp(['Small disturbance, spacing = ' num2str(MV(m)) '; rep # ' num2str(reps)])
        
        % Pick a central reef
        CR = randsample(NumReefs,1);
        [~,I] = sort(Distance(CR,:),'ascend');
        
        % Damage some reefs
        Spacing(1,reps) = MaxSpace;%randi(MaxSpace);
        if MaxSpace == nan
            Damaged = randsample(NumReefs,NumDamaged,0);
        else
            Damaged = I(2:Spacing(1,reps):Spacing(1,reps)*NumDamaged);
        end
        
        % Let them recover
        N = ones(NumSpp,NumReefs);
        for t = 1:1e4
            
            if t == 5
                N(:,Damaged) = 0.4;
            end
            
            for s = 1:NumSpp
                
                % Local-growth driven recovery at a 5% rate
                N(s,:) = min(1,N(s,:)*1.05);
                
                % Background mortality
                Mortrate = ones(1,NumReefs);
                Mortrate(rand(size(Mortrate))<0.4) = 0.85;
                N(s,:) = N(s,:).*Mortrate;
                
                % Quadratic fertilisation success for local production
                LarvalOutput = f.*N(s,:).^2;
                
                % Dispersal-driven recovery
                N(s,:) = min(1,N(s,:) + (LarvalOutput*C(:,:,s)).*(1-N(s,:)));
                SN(s,t) = sum(N(s,:));
            end
            if t > 5 & sum(abs(SN(:,end)-SN(:,1))) < 1
                break
            end
        end
        for s = 1:NumSpp
            Impact(m,s,reps) = sum(abs(SN(s,1) - SN(s,:)));
        end
        clearvars -except reps centroid C Num* Distance CL Impact Spacing MaxSpace k m MV f FS
    end
end
subplot(2,1,2); hold on; d = 0.1;
Q = quantile(Impact,[0.16 0.5 0.84],3)
for m = 1:length(MV)
    for s = 1:NumSpp
        ppp(m) = plot(m+5*s+[-d d nan 0 0 nan -d d],squeeze(Q(m,s,[1 1 1 1 3 3 3 3])),'k','linewidth',2,'color',CL(7-m,:));
        plot(m+5*s,squeeze(Q(m,s,2)),'k.','markersize',12,'color',CL(7-m,:))
    end
end

set(gca,'xtick',[7 12 17 22],'xticklabel',{'Global disperser','Long PLD spawner','Short PLD spawner','Brooder'},'fontsize',FS);
L = legend(ppp,'High correlated disturbances','Low correlation disturbances','No correlation');
set(L,'location','northwest','box','off','fontsize',FS)
title('Small magnitude disturbance','fontsize',FS+2,'fontweight','normal')

load TEMP Impact_Large
Impact_Small = Impact;
save TEMP Impact* MV

Plot_Boxplots
