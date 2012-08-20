function [taus_sorted, P, iterpar] = SubClusterDivision(taus_sorted,P,iterpar)
% SUBCLUSTERDIVISION Divides the two strongest clusters into three subclusters
%
% Author: Mikko Alatossava (CWC/UOULU)
%
% Date: 15.2.2007

LoSConnectionLinks = iterpar.LoSConnectionLinks;
NLoSConnectionLinks = iterpar.NLoSConnectionLinks;
N = iterpar.N;
N_max = max(N(1),N(2));
NumLinks = length(iterpar.UserIndeces);
Scenario = iterpar.Scenario;

%In B1 and B3 NLoS the clusters are not divided in D111
if NLoSConnectionLinks & ~(strcmp(Scenario,'B3')) & ~(strcmp(Scenario,'B1'))            
    tmp_powers = ((P(NLoSConnectionLinks,1:N(2)))./repmat(max((P(NLoSConnectionLinks,1:N(2))),[],2),1,size((P(NLoSConnectionLinks,1:N(2))),2))).';
    
    %find the indeces for the strongest cluster per user
    [MaxValue1, MaxValueIndeces1] = max(tmp_powers,[],1);
    
    %map column indeces MaxValueindeces to real indeces 
    RealIndeces1 = MaxValueIndeces1+[0:size(tmp_powers,1):numel(tmp_powers)-1];
    tmp_powers_zeromax = tmp_powers;
    tmp_powers_zeromax(RealIndeces1) = 0; %set to zero so that it is not the largest anymore
     
    %find the indeces for the second strongest cluster per user
    [MaxValue2, MaxValueIndeces2] = max(tmp_powers_zeromax,[],1);
    RealIndeces2 = MaxValueIndeces2+[0:size(tmp_powers,1):numel(tmp_powers)-1];

    %if the strongest cluster is before the second strongest, the indeces in the second strongest one needs to
    %be shifted by 2. Same applies for vice versa
    Is2BeforeIndex1 = MaxValueIndeces2<MaxValueIndeces1;
    Is1BeforeIndex2 = MaxValueIndeces1<MaxValueIndeces2;

    %when 4 new rows are created before the index in the matrix, the indeces shift.
    %This affects in every row additively, so [0:4:4*size(tmp_powers,2)] is applied
    NewRealIndeces1 = RealIndeces1+Is2BeforeIndex1*2+[0:4:4*size(tmp_powers,2)-1];
    NewRealIndeces2 = RealIndeces2+Is1BeforeIndex2*2+[0:4:4*size(tmp_powers,2)-1];

    %first subcluster has 10 rays, 2nd has 6 and 3rd one 4 rays
    cluster1_subpowers = [10*MaxValue1/20 ;6*MaxValue1/20 ;4*MaxValue1/20]; 
    cluster2_subpowers = [10*MaxValue2/20 ;6*MaxValue2/20 ;4*MaxValue2/20]; 

    %create a matrix to hold all cluster and subclusters and feed the power values into right indeces
    tmp_powers_large = NaN*ones(size(tmp_powers,1)+4,size(tmp_powers,2));
    tmp_powers_large(NewRealIndeces1) = cluster1_subpowers(1,:);
    tmp_powers_large(NewRealIndeces1+1) = cluster1_subpowers(2,:);
    tmp_powers_large(NewRealIndeces1+2) = cluster1_subpowers(3,:);
    tmp_powers_large(NewRealIndeces2) = cluster2_subpowers(1,:);
    tmp_powers_large(NewRealIndeces2+1) = cluster2_subpowers(2,:);
    tmp_powers_large(NewRealIndeces2+2) = cluster2_subpowers(3,:);

    %find out the indeces of the rest of the clusters in tmp_power, basically a setdiff operation
    a=ones(size(tmp_powers));
    b=zeros(size(tmp_powers));
    b(RealIndeces1)=1; b(RealIndeces2)=1;
    OtherIndeces = find(a-b); 

    OtherPowers = tmp_powers(OtherIndeces); %extract those power values to be fed into tmp_powers_large
    OtherPowerLocations = find(isnan(tmp_powers_large)); %indeces in large matrix still needing a value

    tmp_powers_large(OtherPowerLocations) = OtherPowers; %feed the values
    P_nlos = (tmp_powers_large./repmat(sum(tmp_powers_large),size(tmp_powers_large,1),1)).';

    %same for delays
    tmp_taus = (taus_sorted(NLoSConnectionLinks,1:N(2))).';

    %find values corresponding to the max power values
    DelayValue1 = tmp_taus(RealIndeces1);
    DelayValue2 = tmp_taus(RealIndeces2);           

    %insert delay values
    tmp_taus_large = NaN*ones(size(tmp_taus,1)+4,size(tmp_taus,2));
    tmp_taus_large(NewRealIndeces1) = DelayValue1;
    tmp_taus_large(NewRealIndeces1+1) = DelayValue1+5e-9;
    tmp_taus_large(NewRealIndeces1+2) = DelayValue1+10e-9;
    tmp_taus_large(NewRealIndeces2) = DelayValue2;
    tmp_taus_large(NewRealIndeces2+1) = DelayValue2+5e-9;
    tmp_taus_large(NewRealIndeces2+2) = DelayValue2+10e-9;

    OtherDelays = tmp_taus(OtherIndeces);
    tmp_taus_large(OtherPowerLocations) = OtherDelays; %feed the values

    taus_nlos = tmp_taus_large.';

    N(2) = size(taus_nlos,2);
end  

if LoSConnectionLinks
    tmp_powers = ((P(LoSConnectionLinks,1:N(1)))./repmat(max((P(LoSConnectionLinks,1:N(1))),[],2),1,size((P(LoSConnectionLinks,1:N(1))),2))).';
    
    %find the indeces for the strongest cluster per user
    [MaxValue1, MaxValueIndeces1] = max(tmp_powers,[],1);
    
    %map column indeces MaxValueindeces to real indeces 
    RealIndeces1 = MaxValueIndeces1+[0:size(tmp_powers,1):numel(tmp_powers)-1];
    tmp_powers_zeromax = tmp_powers;
    tmp_powers_zeromax(RealIndeces1) = 0; %set to zero so that it is not the largest anymore
     
    %find the indeces for the second strongest cluster per user
    [MaxValue2, MaxValueIndeces2] = max(tmp_powers_zeromax,[],1);
    RealIndeces2 = MaxValueIndeces2+[0:size(tmp_powers,1):numel(tmp_powers)-1];

    %if the strongest cluster is before the second strongest, the indeces in the second strongest one needs to
    %be shifted by 2. Same applies for vice versa
    Is2BeforeIndex1 = MaxValueIndeces2<MaxValueIndeces1;
    Is1BeforeIndex2 = MaxValueIndeces1<MaxValueIndeces2;

    %when 4 new rows are created before the index in the matrix, the indeces shift.
    %This affects in every row additively, so [0:4:4*size(tmp_powers,2)] is applied
    NewRealIndeces1 = RealIndeces1+Is2BeforeIndex1*2+[0:4:4*size(tmp_powers,2)-1];
    NewRealIndeces2 = RealIndeces2+Is1BeforeIndex2*2+[0:4:4*size(tmp_powers,2)-1];

    %first subcluster has 10 rays, 2nd has 6 and 3rd one 4 rays
    cluster1_subpowers = [10*MaxValue1/20 ;6*MaxValue1/20 ;4*MaxValue1/20]; 
    cluster2_subpowers = [10*MaxValue2/20 ;6*MaxValue2/20 ;4*MaxValue2/20]; 

    %create a matrix to hold all cluster and subclusters and feed the power values into right indeces
    tmp_powers_large = NaN*ones(size(tmp_powers,1)+4,size(tmp_powers,2));
    tmp_powers_large(NewRealIndeces1) = cluster1_subpowers(1,:);
    tmp_powers_large(NewRealIndeces1+1) = cluster1_subpowers(2,:);
    tmp_powers_large(NewRealIndeces1+2) = cluster1_subpowers(3,:);
    tmp_powers_large(NewRealIndeces2) = cluster2_subpowers(1,:);
    tmp_powers_large(NewRealIndeces2+1) = cluster2_subpowers(2,:);
    tmp_powers_large(NewRealIndeces2+2) = cluster2_subpowers(3,:);

    %find out the indeces of the rest of the clusters in tmp_power, basically a setdiff operation
    a=ones(size(tmp_powers));
    b=zeros(size(tmp_powers));
    b(RealIndeces1)=1; b(RealIndeces2)=1;
    OtherIndeces = find(a-b); 

    OtherPowers = tmp_powers(OtherIndeces); %extract those power values to be fed into tmp_powers_large
    OtherPowerLocations = find(isnan(tmp_powers_large)); %indeces in large matrix still needing a value

    tmp_powers_large(OtherPowerLocations) = OtherPowers; %feed the values
    P_los = (tmp_powers_large./repmat(sum(tmp_powers_large),size(tmp_powers_large,1),1)).';
    
    
    %same for delays
    tmp_taus = (taus_sorted(LoSConnectionLinks,1:N(1))).';

    %find values corresponding to the max power values
    DelayValue1 = tmp_taus(RealIndeces1);
    DelayValue2 = tmp_taus(RealIndeces2);           

    %insert delay values
    tmp_taus_large = NaN*ones(size(tmp_taus,1)+4,size(tmp_taus,2));
    tmp_taus_large(NewRealIndeces1) = DelayValue1;
    tmp_taus_large(NewRealIndeces1+1) = DelayValue1+5e-9;
    tmp_taus_large(NewRealIndeces1+2) = DelayValue1+10e-9;
    tmp_taus_large(NewRealIndeces2) = DelayValue2;
    tmp_taus_large(NewRealIndeces2+1) = DelayValue2+5e-9;
    tmp_taus_large(NewRealIndeces2+2) = DelayValue2+10e-9;

    OtherDelays = tmp_taus(OtherIndeces);
    tmp_taus_large(OtherPowerLocations) = OtherDelays; %feed the values

    taus_los = tmp_taus_large.';

    N(1) = size(taus_los,2);
end

N_max = max(N(1),N(2));
taus_sorted_final = NaN*ones(NumLinks,N_max);
P_final = NaN*ones(NumLinks,N_max);

if LoSConnectionLinks
    taus_sorted_final(LoSConnectionLinks,1:N(1)) = taus_los;
    P_final(LoSConnectionLinks,1:N(1)) = P_los;
else
    taus_sorted_final(LoSConnectionLinks,1:N(1)) = taus_sorted(LoSConnectionLinks,1:N(1));
    P_final(LoSConnectionLinks,1:N(1)) = P(LoSConnectionLinks,1:N(1));
end

if NLoSConnectionLinks & ~(strcmp(Scenario,'B3')) & ~(strcmp(Scenario,'B1'));
    taus_sorted_final(NLoSConnectionLinks,1:N(2)) = taus_nlos;
    P_final(NLoSConnectionLinks,1:N(2)) = P_nlos;
else
    taus_sorted_final(NLoSConnectionLinks,1:N(2)) = taus_sorted(NLoSConnectionLinks,1:N(2));
    P_final(NLoSConnectionLinks,1:N(2)) = P(NLoSConnectionLinks,1:N(2));
end

taus_sorted = taus_sorted_final;
P = P_final;
iterpar.N = N;