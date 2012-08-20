function bulk_parameters = generate_bulk_par(wimpar,linkpar,antpar,fixpar)
%GENERATE_BULK_PAR Generation of WIM2 interim channel model parameters
%   [BULK_PAR]=GENERATE_BULK_PAR(WIMPAR,LINKPAR,ANTPAR,FIXPAR) generates the
%   "bulk" parameters according to WINNER D5.4 with some Phase II modifications.
%   For explanation of the input structs, see WIMPARSET, LINKPARSET, and ANTPARSET.
%   Denoting with K the number of links, N the number of paths,
%   M the number of subpaths, the fields BULK_PAR are as follows:
%
%   delays           - path delays in seconds [KxN]
%   path_powers      - relative path powers [KxN]
%   aods             - angles of departure in degrees over (-180,180) [KxNxM]
%   aoas             - angles of arrival in degrees over (-180,180) [KxNxM]
%   subpath_phases   - random phases for subpaths in degrees over (0,360) [KxNxM]
%   path_losses      - path losses in linear scale [Kx1]
%   MsBsDistance     - distances between MSs and BSs in meters [1xK]
%   shadow_fading    - shadow fading losses in linear scale [Kx1]
%   propag_condition -whether the user is in LoS condition (1) or in nlos (0)
%   sigmas           -correlation coefficients fo large scale parameters
%
%   In addition, when users with LoS condition exists (in addition to the above):
%   Kcluster        - K factors for all links [Kx1]
%   Phi_LOS         - random phases for LOS paths in degrees over (-180,180) [Kx1]
%
%   In addition, when users wimpar.PolarisedArrays is 'yes' (in addition to the above):
%   xprV            -vertical xpr values, [KxNxM]
%   xprH            -horizontal xpr values, [KxNxM]
%
%   In addition, when users in B5 scenario exist (in addition to the above):
%   scatterer_freq  -Doppler frequency for scatterers, [KxNxM]
%
%   Ref. [1]: D1.1.1 V1.0, "WINNER II interim channel models"
%        [2]: 3GPP TR 25.996 v6.1.0 (2003-09)
%        [3]: D. Reed et. al, "Spatial Channel Models for Multi-antenna
%             Systems, ..., 2003.
%
%   See also WIM.

%   Authors of model versions:
%
%   WINNER Phase II interim (WIM2i): Pekka Ky�sti (EBIT), Lassi Hentil� (EBIT),
%   Marko Milojevic (TUI), Mikko Alatossava (CWC/UOULU)
%
%   WINNER Phase I (WIM): Daniela Laselva (EBIT), Marko Milojevic (TUI),
%   Pekka Ky�sti (EBIT), Lassi Hentil� (EBIT)
%
%   SCM/SCME: Jari Salo (HUT), Daniela Laselva (EBIT), Giovanni Del Galdo (TUI),
%   Marko Milojevic (TUI), Pekka Ky�sti (EBIT), Christian Schneider (TUI),
%   Zhiwen Wu(BUPT),Yu Zhang(BUPT),Jianhua Zhang(BUPT),Guangyi Liu(CMCC)

% Prevent using CDL models (to be removed when CDL model tables are agreed)
%if strcmpi(wimpar.FixedPdpUsed,'yes') | strcmpi(wimpar.FixedAnglesUsed,'yes')
%error('CDL models not supported. Set wimpar FixedPdpUsed=NO and FixedAnglesUsed=NO')
%end
%   Update InH AOA/AOD distribution to Laplacian      25.9 2008 Zhiwen WU

% Number of scenarios in WINNER II channel models
NumOfScenarios = 6;    % should be equal to the number in ScenarioMapping.m

% extract certain parameters from the input structs
ScenarioVector        = linkpar.ScenarioVector;
MsBsDistance          = linkpar.MsBsDistance;

bulk_parameters = struct_generation(1, 1, wimpar, linkpar,1 , 'Initialization');

for ScenIndex = 1:NumOfScenarios
    %save parameters that vary between different for loop iterations to iterpar
    iterpar.Scenario = ScenarioMapping(ScenIndex); %map scenario from numerical value to letters
    iterpar.UserIndeces = find(ScenarioVector==ScenIndex); %links (user) that are in a certain scenario correspongind to ScenIndex

    if iterpar.UserIndeces
        switch iterpar.Scenario
            % Geometric based stochastic models
            case {'A2', 'B1', 'B4','C1', 'C2', 'D1'}
                bulk_parameters_iter = stochastic(wimpar,linkpar,antpar,fixpar,iterpar);

                % B5, Static feeder scenario
            case {'B5a','B5b','B5c','B5f'}
                bulk_parameters_iter = static(wimpar,linkpar,antpar,fixpar,iterpar);

        end     % end of user parameter generation main program

        bulk_parameters = struct_generation(bulk_parameters, bulk_parameters_iter, wimpar, linkpar, iterpar,'Iteration');
        clear bulk_parameters_iter iterpar
    end

end

bulk_parameters = struct_generation(bulk_parameters, 1, wimpar, linkpar, 1, 'Refinement');

% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to generate bulk parameters
% See [1, Sec.3.1.7].
function bulk_parameters = stochastic(wimpar,linkpar,antpar,fixpar,iterpar)

%-- STEP 1 --%
% extract certain parameters from the input structs
PolarisedArrays       = wimpar.PolarisedArrays;
M                     = wimpar.NumSubPathsPerPath;
DelaySamplingInterval = wimpar.DelaySamplingInterval;
FixedPdpUsed          = wimpar.FixedPdpUsed;

%-- STEP 2 --%
% extract the number of users from the first field of linkpar struct
UserIndeces = iterpar.UserIndeces;
MsBsDistance = linkpar.MsBsDistance(UserIndeces);
NumLinks     = length(UserIndeces);
Scenario     = iterpar.Scenario;

switch Scenario
    case {'A2', 'B1', 'C1', 'C2', 'D1',}
        evalstr = sprintf('iterpar.LoS = fixpar.%s.LoS;',Scenario);
        eval(evalstr);
        evalstr = sprintf('iterpar.NLoS = fixpar.%s.NLoS;',Scenario);
        eval(evalstr);

        %number of cluster differs in LoS and NLoS situations
        N(1) = iterpar.LoS.NumClusters;
        N(2) = iterpar.NLoS.NumClusters;

    case {'B4'}
        evalstr = sprintf('iterpar.NLoS = fixpar.%s.NLoS;',Scenario);
        eval(evalstr);
        N(1) = 0;
        N(2) = iterpar.NLoS.NumClusters;
end
N_max = max(N(1),N(2));
iterpar.N = N;

% check that M = 20
if (M ~= 20)
    M=20;
    warning('MATLAB:NumSubPathsChanged','NumSubPathsPerPath is not 20! Using NumSubPathsPerPath=20 instead.')
end

% make sure that user-specific parameters are row vectors
ThetaBs      = linkpar.ThetaBs(UserIndeces);
ThetaMs      = linkpar.ThetaMs(UserIndeces);

%extract line of sight probability
if strcmp(wimpar.UseManualPropCondition,'yes')
    PropagCondition = linkpar.PropagConditionVector(UserIndeces);
    switch Scenario
        case {'C4'} %all nlos
            PropagCondition = zeros(1,NumLinks);
    end
else
    switch Scenario
        case {'A2', 'B1', 'C1', 'C2', 'D1'}
            PropagCondition = LOSprobability(wimpar,linkpar,fixpar,iterpar);
        case {'B4'}
            PropagCondition = zeros(1,NumLinks);
    end
end

%indeces of LoS/NLoS links and the amount of them
iterpar.PropagCondition = PropagCondition;
LoSConnectionLinks = find(PropagCondition); iterpar.LoSConnectionLinks = LoSConnectionLinks;
NumLoSConnectionLinks = length(LoSConnectionLinks); iterpar.NumLoSConnectionLinks = NumLoSConnectionLinks;
NLoSConnectionLinks = find(PropagCondition==0); iterpar.NLoSConnectionLinks = NLoSConnectionLinks;
NumNLoSConnectionLinks = length(NLoSConnectionLinks); iterpar.NumNLoSConnectionLinks = NumNLoSConnectionLinks;

if NumLoSConnectionLinks == 0
    N(1) = 0;
end
if NumNLoSConnectionLinks == 0
    N(2) = 0;
end
N_max = max(N(1),N(2));
iterpar.N = N;


%FS, STEP 1-2, all users exhibit bad urban effect (long delays)
if strcmp(Scenario,'B2') | strcmp(Scenario,'C3')
    if strcmp(Scenario,'B2') 
        MsScatBsDist = sort(1000-(1000-300)*rand(NumNLoSConnectionLinks,2),2); %generate 2 scatterers for each user with distances from [1, table 4-3]
        FSLoss = 4; %power loss due to excess delay dB per us
    else
        MsScatBsDist = sort(3000-(3000-600)*rand(NumNLoSConnectionLinks,2),2); %generate 2 scatterers for each user with distances from [1, table 4-3]
        FSLoss = 2; %power loss due to excess delay dB per us
    end

    NumFSConnectionLinks = NumNLoSConnectionLinks;
    FSConnectionLinks = NLoSConnectionLinks;
    NumFSPaths = 2; %two last clusters for each path are created as FS (Far Scatter) clusters
else
    FSConnectionLinks = [];
    FSPaths = 0;
end



%-- STEP 3 --%
% employ the user-defined path loss model
%if isequal(lower(wimpar.PathLossModelUsed),'yes')
    [path_losses, linkpar, fixpar, iterpar] = feval(wimpar.PathLossModel,wimpar,linkpar,fixpar,iterpar);
    path_losses = 10.^(-path_losses(:)/10);    % a (NumLinks x 1) vector
%else
%    path_losses=NaN*ones(1,length(iterpar.UserIndeces));
%end

%-- STEP 4 --%
% Generation of correlated DS, AS's and SF for all users
% This step takes into account channel scenario automatically
sigmas = LScorrelation(wimpar,linkpar,fixpar,iterpar);
sigma_asD = sigmas(:,1);
sigma_asA = sigmas(:,2);
sigma_ds  = sigmas(:,3);
sigma_sf  = sigmas(:,4);
sigma_kf  = sigmas(:,5);

% generate vehicle penetration shadowing
sigma_sf_o2v = 10.^(0.1*5*randn(1,NumLinks)); %log-normal dB

%-- STEP 5 --%
%% Generate delays in a (NumLinks x N) matrix %%
% The unit of taus is seconds
if strcmpi(FixedPdpUsed,'no')
    sigma_ds = repmat(sigma_ds,1,N_max);             % delay spreads for all clusers/users
    taus = NaN*ones(NumLinks,N_max);
    taus_sorted = NaN*ones(NumLinks,N_max);
    taus_los = NaN*ones(NumLinks,N_max);
    taus_rounded = NaN*ones(NumLinks,N_max);

%     switch upper(Scenario)      % See distributions in [1, table 4-5]
% 
%         case {'A1','A2','B3','B4','C1','C2','C3','C4','D1','D2A'}
            if LoSConnectionLinks
                taus(LoSConnectionLinks,1:N(1)) = sort(-iterpar.LoS.r_DS*sigma_ds(LoSConnectionLinks,1:N(1)).*log(rand(NumLoSConnectionLinks,N(1))),2);  % Exp [1, eq. 4.1]
            end

            if NLoSConnectionLinks
                taus(NLoSConnectionLinks,1:N(2)) = sort(-iterpar.NLoS.r_DS*sigma_ds(NLoSConnectionLinks,1:N(2)).*log(rand(NumNLoSConnectionLinks,N(2))),2);  % Exp [1, eq. 4.1]
            end


%     end     % end switch

    taus_sorted(LoSConnectionLinks,1:N(1))  = taus(LoSConnectionLinks,1:N(1)) - repmat(taus(LoSConnectionLinks,1),1,N(1));       % normalize min. delay to zero
    taus_sorted(NLoSConnectionLinks,1:N(2))  = taus(NLoSConnectionLinks,1:N(2)) - repmat(taus(NLoSConnectionLinks,1),1,N(2));

    %FS STEP 4
    if FSConnectionLinks
        taus_sorted(FSConnectionLinks,N(2)-NumFSPaths+1:N(2)) = 0;
        ExcessDelayLoss = ((MsScatBsDist - repmat(MsBsDistance(FSConnectionLinks).',1,NumFSPaths))./ ...
            (repmat(3e8,NumFSConnectionLinks,NumFSPaths))).*1e6*FSLoss;
    end

    %in case of los, extra factor. Not be used in clusterpower calculation
    taus_without_los_factor = taus_sorted; %need to extract taus_sorted to taus_los, since the next step is not applied for the taus that is given as input parameter to powers generation
    if FSConnectionLinks
        taus_sorted(FSConnectionLinks,end-1:end) = MsScatBsDist./3e8; 
    end
    if NumLoSConnectionLinks
        K_factors = sigma_kf';
        K_factors_dB = 10*log10(abs(sigma_kf))';
        ConstantD = 0.7705-0.0433.*K_factors_dB(LoSConnectionLinks) + ...
            0.0002.*K_factors_dB(LoSConnectionLinks).^2 + 0.000017.*K_factors_dB(LoSConnectionLinks).^3; %[1, eq.4.3]
        taus_sorted(LoSConnectionLinks,1:N(1)) = taus_sorted(LoSConnectionLinks,1:N(1))./repmat(ConstantD.',1,N(1));
    end

else    % use fixed delays from a table
    [taus_los,Pprime_los,Kcluster_los,...
        taus_nlos,Pprime_nlos,Kcluster_nlos] = fixedPdp(wimpar,iterpar);    % the same for each link

    taus_sorted = NaN*ones(NumLinks,N_max);

    taus_sorted(LoSConnectionLinks,1:N(1)) = repmat(taus_los,NumLoSConnectionLinks,1);
    taus_sorted(NLoSConnectionLinks,1:N(2)) = repmat(taus_nlos,NumNLoSConnectionLinks,1);
    K_factors_dB(LoSConnectionLinks) = Kcluster_los(1,:); %dB
    K_factors(LoSConnectionLinks) = 10.^(K_factors_dB(LoSConnectionLinks)/10); %linear

end

% Rounding to delay grid
if (DelaySamplingInterval>0)
    taus_rounded = DelaySamplingInterval*floor(taus_sorted/DelaySamplingInterval + 0.5);
else
    taus_rounded = taus_sorted;
end

% end of delay generation

%-- STEP 6 --%
%% Determine random average powers in a (NumLinks x N) matrix %%
if strcmpi(FixedPdpUsed,'no')
    if LoSConnectionLinks
        ksi_LoS = randn(NumLoSConnectionLinks,N(1))*iterpar.LoS.LNS_ksi;           % per-path shadowing
    end
    if NLoSConnectionLinks
        ksi_NLoS = randn(NumNLoSConnectionLinks,N(2))*iterpar.NLoS.LNS_ksi;           % per-path shadowing
    end
    P = NaN*ones(NumLinks,N_max);
    Pprime = NaN*ones(NumLinks,N_max);

    % See distributions in [1, table 4-5]
    %for LoS links, with exponential distribution of delays [1, eq 4.3]
    if LoSConnectionLinks
        Pprime(LoSConnectionLinks,1:N(1)) =  exp(-taus_without_los_factor(LoSConnectionLinks,1:N(1)).*((iterpar.LoS.r_DS-1)./(iterpar.LoS.r_DS.*sigma_ds(LoSConnectionLinks,1:N(1))))).*10.^(-ksi_LoS/10);
    end
    if NLoSConnectionLinks
        %for NLoS links, with exponential distribution of delays [1, eq 4.3]
        Pprime(NLoSConnectionLinks,1:N(2)) =  exp(-taus_without_los_factor(NLoSConnectionLinks,1:N(2)).*((iterpar.NLoS.r_DS-1)./(iterpar.NLoS.r_DS.*sigma_ds(NLoSConnectionLinks,1:N(2))))).*10.^(-ksi_NLoS/10);

        %FS STEP 5
        if FSConnectionLinks
            Pprime(FSConnectionLinks,N(2)+1-NumFSPaths:N(2)) = Pprime(FSConnectionLinks,N(2)+1-NumFSPaths:N(2)) ...
                .* 10.^(-ExcessDelayLoss./10);
        end
    end

    %for LoS links
    if NumLoSConnectionLinks
        %temporary P_tmp that is used to replace P(LoSConnectionLinks,1:N(1)) after the angular directions have been created
        P_tmp = Pprime(LoSConnectionLinks,1:N(1))./repmat(sum(Pprime(LoSConnectionLinks,1:N(1)),2),1,N(1));
        %Kfactor calculations are here only for angular domain use
        SpecularRayPower = K_factors(LoSConnectionLinks)./(K_factors(LoSConnectionLinks)+1); %[1, eq.4.8]
        DiracVector = zeros(NumLoSConnectionLinks,N(1));
        DiracVector(:,1)=1;
        P(LoSConnectionLinks,1:N(1)) = repmat(1./(1+K_factors(LoSConnectionLinks).'),1,N(1)).*(Pprime(LoSConnectionLinks,1:N(1)) ...
            ./repmat(sum(Pprime(LoSConnectionLinks,1:N(1)),2),1,N(1))) + ...
            DiracVector.*repmat(SpecularRayPower.',1,N(1)); %[1, eq.4.9]
    end

    %for NLoS links
    if NumNLoSConnectionLinks
        P(NLoSConnectionLinks,1:N(2)) = Pprime(NLoSConnectionLinks,1:N(2))./repmat(sum(Pprime(NLoSConnectionLinks,1:N(2)),2),1,N(2));
    end

else    % use fixed powers from a table
    [taus_los,Pprime_los,Kcluster_los,...
        taus_nlos,Pprime_nlos,Kcluster_nlos] = fixedPdp(wimpar,iterpar);    % the same for each link

    % Replace number of paths by number of tabulated paths
    P = NaN*ones(NumLinks,N_max);
    P(LoSConnectionLinks,1:N(1)) = repmat(Pprime_los,NumLoSConnectionLinks,1);
    P(NLoSConnectionLinks,1:N(2)) = repmat(Pprime_nlos,NumNLoSConnectionLinks,1);
    P(LoSConnectionLinks,1:N(1)) = P(LoSConnectionLinks,1:N(1))./repmat(sum(P(LoSConnectionLinks,1:N(1)),2),1,N(1));
    P(NLoSConnectionLinks,1:N(2)) = P(NLoSConnectionLinks,1:N(2))./repmat(sum(P(NLoSConnectionLinks,1:N(2)),2),1,N(2));
end

%-- STEP 7 --%
%% Determine AoDs / AoAs %%
offset = [0.0447 0.1413 0.2492 0.3715 0.5129 0.6797 0.8844 1.1481 1.5195 2.1551]; % [1, Table 4-1] +/- offset angles, resulting Laplacian APS, with rms AS = 1 deg

if strcmpi(wimpar.FixedAnglesUsed,'no')

    [offset_matrix_AoD, offset_matrix_AoA] = offset_matrix_generation(offset,iterpar);

    AoDPrimer = NaN*ones(NumLinks,N_max);
    AoD_path = NaN*ones(NumLinks,N_max);
    AoAPrimer = NaN*ones(NumLinks,N_max);
    AoA_path = NaN*ones(NumLinks,N_max);

    %pick a correct scaling factor C in equation [1, Eq. 4.10]
    ScalingFactorC_matrix = NaN*ones(NumLinks,N_max);
    % Table of constant C in [1, step 7]
    switch upper(Scenario)
        case{'B1','C1','C2','B4','D1'}
            ConstantC = [4 5 8 10 11 12 14 15 16 19 20;...
            0.779 0.860 1.018 1.090 1.123 1.146 1.190 1.211 1.226 1.273 1.289];
        case{'A2'}
            ConstantC = [4 5 8 10 11 12 14 15 16 19 20;...
            0.779 0.860 1.018 1.090 1.123 1.146 1.190 1.434 1.226 1.501 1.289];
    end;
    if LoSConnectionLinks
        ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)) = ConstantC(2,find(ConstantC(1,:)==N(1)));
        %ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)) = ConstantCMapping(N(1));
        ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)) = ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)).*repmat((1.1035-0.028.*K_factors_dB(LoSConnectionLinks).'- ...
            0.002.*(K_factors_dB(LoSConnectionLinks).').^2+0.0001.*(K_factors_dB(LoSConnectionLinks).').^3),1,N(1));  %[1, eq.4.11]

        switch upper(Scenario)
            case{'B1','C1','C2','B4','D1'}
                ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)) = ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)).*repmat((1.1035-0.028.*K_factors_dB(LoSConnectionLinks).'- ...
                0.002.*(K_factors_dB(LoSConnectionLinks).').^2+0.0001.*(K_factors_dB(LoSConnectionLinks).').^3),1,N(1));  %[1, eq.4.11]
                AoDPrimer(LoSConnectionLinks,1:N(1)) = (2*repmat(sigma_asD(LoSConnectionLinks),1,N(1))/1.4.*sqrt(-log(P(LoSConnectionLinks,1:N(1))./repmat(max(P(LoSConnectionLinks,1:N(1)),[],2),1,N(1)))))./ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)); %[1, eq. 4.10]
                AoAPrimer(LoSConnectionLinks,1:N(1)) = (2*repmat(sigma_asA(LoSConnectionLinks),1,N(1))/1.4.*sqrt(-log(P(LoSConnectionLinks,1:N(1))./repmat(max(P(LoSConnectionLinks,1:N(1)),[],2),1,N(1)))))./ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)); %[1, eq. 4.10]
            case{'A2'}
                ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)) = ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)).*repmat((0.9275+0.0439.*K_factors_dB(LoSConnectionLinks).'- ...
                0.0071.*(K_factors_dB(LoSConnectionLinks).').^2+0.0002.*(K_factors_dB(LoSConnectionLinks).').^3),1,N(1));  %[1, eq.4.11]
                AoDPrimer(LoSConnectionLinks,1:N(1)) = -(repmat(sigma_asD(LoSConnectionLinks),1,N(1)).*(log(P(LoSConnectionLinks,1:N(1))./repmat(max(P(LoSConnectionLinks,1:N(1)),[],2),1,N(1)))))./ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)); %[1, eq. 4.10]
                AoAPrimer(LoSConnectionLinks,1:N(1)) = -(repmat(sigma_asA(LoSConnectionLinks),1,N(1)).*(log(P(LoSConnectionLinks,1:N(1))./repmat(max(P(LoSConnectionLinks,1:N(1)),[],2),1,N(1)))))./ScalingFactorC_matrix(LoSConnectionLinks,1:N(1)); %[1, eq. 4.10]

        end;
        
        X_AoD = (2*round(rand(NumLoSConnectionLinks,N(1)))-1);
        Y_AoD = repmat(sigma_asD(LoSConnectionLinks),1,N(1))/1.4/5.*randn(NumLoSConnectionLinks,N(1));
        AoD_path(LoSConnectionLinks,1:N(1)) = (AoDPrimer(LoSConnectionLinks,1:N(1)).*X_AoD + Y_AoD - ...
            (repmat(AoDPrimer(LoSConnectionLinks,1),1,N(1))*X_AoD(1) + Y_AoD(1) - repmat(ThetaBs(LoSConnectionLinks).',1,N(1)))); %[1, eq. 4.13]


        X_AoA = (2*round(rand(NumLoSConnectionLinks,N(1)))-1);
        Y_AoA = repmat(sigma_asA(LoSConnectionLinks),1,N(1))/1.4/5.*randn(NumLoSConnectionLinks,N(1));
        AoA_path(LoSConnectionLinks,1:N(1)) = (AoAPrimer(LoSConnectionLinks,1:N(1)).*X_AoA + Y_AoA - ...
            (repmat(AoAPrimer(LoSConnectionLinks,1),1,N(1))*X_AoA(1) + Y_AoA(1) - repmat(ThetaBs(LoSConnectionLinks).',1,N(1)))); %[1, eq. 4.13]


    end

    if NLoSConnectionLinks
        ScalingFactorC_matrix(NLoSConnectionLinks,1:N(2)) = ConstantC(2,find(ConstantC(1,:)==N(2)));
        %ScalingFactorC_matrix(NLoSConnectionLinks,1:N(2)) = ConstantCMapping(N(2));
         switch upper(Scenario)
            case{'B1','C1','C2','B4','D1'}
                AoDPrimer(NLoSConnectionLinks,1:N(2)) = (2*repmat(sigma_asD(NLoSConnectionLinks),1,N(2))/1.4.*sqrt(-log(P(NLoSConnectionLinks,1:N(2))./repmat(max(P(NLoSConnectionLinks,1:N(2)),[],2),1,N(2)))))./ScalingFactorC_matrix(NLoSConnectionLinks,1:N(2)); %[1, eq. 4.10]
                AoAPrimer(NLoSConnectionLinks,1:N(2)) = (2*repmat(sigma_asA(NLoSConnectionLinks),1,N(2))/1.4.*sqrt(-log(P(NLoSConnectionLinks,1:N(2))./repmat(max(P(NLoSConnectionLinks,1:N(2)),[],2),1,N(2)))))./ScalingFactorC_matrix(NLoSConnectionLinks,1:N(2)); %[1, eq. 4.10]
            case{'A2'}
                AoDPrimer(NLoSConnectionLinks,1:N(2)) = -(repmat(sigma_asD(NLoSConnectionLinks),1,N(2)).*(log(P(NLoSConnectionLinks,1:N(2))./repmat(max(P(NLoSConnectionLinks,1:N(2)),[],2),1,N(2)))))./ScalingFactorC_matrix(NLoSConnectionLinks,1:N(2)); %[1, eq. 4.10]
                AoAPrimer(NLoSConnectionLinks,1:N(2)) = -(repmat(sigma_asA(NLoSConnectionLinks),1,N(2)).*(log(P(NLoSConnectionLinks,1:N(2))./repmat(max(P(NLoSConnectionLinks,1:N(2)),[],2),1,N(2)))))./ScalingFactorC_matrix(NLoSConnectionLinks,1:N(2)); %[1, eq. 4.10]
        end;
        
        AoD_path(NLoSConnectionLinks,1:N(2)) = AoDPrimer(NLoSConnectionLinks,1:N(2)).* (2*round(rand(NumNLoSConnectionLinks,N(2)))-1) + ...
            repmat(sigma_asD(NLoSConnectionLinks),1,N(2))/1.4/5.*randn(NumNLoSConnectionLinks,N(2)) + ...
            repmat(ThetaBs(NLoSConnectionLinks).',1,N(2)); %[1, eq. 4.12]

        AoA_path(NLoSConnectionLinks,1:N(2)) = AoAPrimer(NLoSConnectionLinks,1:N(2)).* (2*round(rand(NumNLoSConnectionLinks,N(2)))-1) + ...
            repmat(sigma_asA(NLoSConnectionLinks),1,N(2))/1.4/5.*randn(NumNLoSConnectionLinks,N(2)) + ...
            repmat(ThetaMs(NLoSConnectionLinks).',1,N(2)); %[1, eq. 4.12]

    end

    AoD_tmp = repmat(reshape(AoD_path.',1,N_max*NumLinks),M,1);
    theta_nm_aod = AoD_tmp + offset_matrix_AoD; %[1, eq. 4.14] M x(NxNumLinks) matrix

    AoA_tmp = repmat(reshape(AoA_path.',1,N_max*NumLinks),M,1);
    theta_nm_aoa = AoA_tmp + offset_matrix_AoA; %[1, eq. 4.14] M x(NxNumLinks) matrix


    %-- STEP 8 --%
    % Pair AoA rays randomly with AoD rays (within a cluster)
    [dummy h]           = sort(rand(M,N_max*NumLinks),1);       % create N*NumLinks random permutations of integers [1:M]
    inds                = h+repmat([1:M:M*N_max*NumLinks],M,1)-1;
    theta_nm_aoa        = theta_nm_aoa(inds);    % random permutation of columns, a (M x N*NumLinks) matrix


else % use fixed AoD/AoAs (without random pairing of subpaths)

    % Determine AoDs %%
    [AoD_path_los,iterpar.PerClusterAS_D,...
        AoD_path_nlos,iterpar.NLoS.PerClusterAS_D] = fixedAods(wimpar,iterpar);    % the same for each link

    %Determine AoAs
    [AoA_path_los,iterpar.PerClusterAS_A,...
        AoA_path_nlos,iterpar.NLoS.PerClusterAS_A] = fixedAoas(wimpar,iterpar);    % the same for each link

    AoD_path = NaN*ones(NumLinks,N_max);
    AoA_path = NaN*ones(NumLinks,N_max);

    if LoSConnectionLinks
        AoD_path(LoSConnectionLinks,1:N(1)) = repmat(AoD_path_los,NumLoSConnectionLinks,1);
        AoA_path(LoSConnectionLinks,1:N(1)) = repmat(AoA_path_los,NumLoSConnectionLinks,1);
    end

    if NLoSConnectionLinks
        AoD_path(NLoSConnectionLinks,1:N(2)) = repmat(AoD_path_nlos,NumNLoSConnectionLinks,1);
        AoA_path(NLoSConnectionLinks,1:N(2)) = repmat(AoA_path_nlos,NumNLoSConnectionLinks,1);
    end

    AoD_tmp = repmat(reshape(AoD_path.',1,N_max*NumLinks),M,1);
    AoA_tmp = repmat(reshape(AoA_path.',1,N_max*NumLinks),M,1);

    [offset_matrix_AoD, offset_matrix_AoA] = offset_matrix_generation(offset,iterpar);

    %apply offset matrix
    % NOTE! now array orientation parameter ThetaBs is disabled and AoD is always
    % like in CDL model tables, 17.5.2006 PekKy.
    theta_nm_aod    = AoD_tmp + offset_matrix_AoD;      % a (M x (NumLinks*N)) matrix
    theta_nm_aoa    = AoA_tmp + offset_matrix_AoA;      % a (M x (NumLinks*N)) matrix
end


% Values of theta_nm_aoa and theta_nm_aod may be outside (-180,180).
% Wrapping of angles to range (-180,180)
theta_nm_aoa = prin_value(theta_nm_aoa);
theta_nm_aod = prin_value(theta_nm_aod);


% put AoDs and AoAs into a 3D-array with dims [NumLinks N M]
theta_nm_aod=reshape(theta_nm_aod,M,N_max,NumLinks);
theta_nm_aod=permute(theta_nm_aod,[3 2 1]);
theta_nm_aoa=reshape(theta_nm_aoa,M,N_max,NumLinks);
theta_nm_aoa=permute(theta_nm_aoa,[3 2 1]);

%-- STEP 10a --%
phi = 360*rand(NumLinks,N_max,M);        % random phases for all users, Uni(0,360)

%set to NaN those that are not valid
if N(1) < N(2)
    phi(LoSConnectionLinks,end+1-(N(2)-N(1)):end,:) = NaN;
elseif N(1) < N(2)
    phi(NLoSConnectionLinks,end+1-(N(1)-N(2)):end,:) = NaN;
end

%replace the kfactor related powers with powers independent of Kfactor. Kfactor information will be used in wim_core-function
if NumLoSConnectionLinks & strcmp(FixedPdpUsed,'no')
    P(LoSConnectionLinks,1:N(1)) = P_tmp; %P_tmp exists only in LoS case
end

%%% Output generation %%%
bulk_parameters=struct( 'delays',taus_rounded,...
    'path_powers',P,...
    'aods',theta_nm_aod,...         % in degrees
    'aoas',theta_nm_aoa,...         % in degrees
    'subpath_phases',phi,...        % in degrees
    'path_losses',path_losses,...   % in linear scale
    'MsBsDistance',MsBsDistance,... % This output is needed since the originally generated MsBsDistances are fitted inside the applicabity ranges of the Scenarios
    'shadow_fading',sigma_sf,...    % in linear scale
    'sigmas',sigmas, ...
    'propag_condition', PropagCondition.',...
    'LoS02ILinks', linkpar.LoS02ILinks,...
    'LoS02VLinks', linkpar.LoS02VLinks,...
    'NLoS02ILinks', linkpar.NLoS02ILinks,...
    'NLoS02VLinks', linkpar.NLoS02VLinks,...
    'o2v_shadow_fading',sigma_sf_o2v,...
    'user_indeces',UserIndeces.',...
    'Kcluster',sigmas(:,5)');

if strcmpi(FixedPdpUsed,'no')
    Phi_LOS = NaN*ones(NumLinks,1);
    Phi_LOS(LoSConnectionLinks) = 360*(rand(NumLoSConnectionLinks,1)-0.5);
    bulk_parameters.Phi_LOS = Phi_LOS;

elseif strcmpi(FixedPdpUsed,'yes')
    % set the LOS phase randomly
    Phi_LOS = NaN*ones(NumLinks,1);
    Phi_LOS(LoSConnectionLinks) = 360*(rand(NumLoSConnectionLinks,size(K_factors,1))-0.5);

    bulk_parameters.Phi_LOS = Phi_LOS;
end

%-- STEP 9 - STEP10b --%
if strcmpi(PolarisedArrays,'yes')

    % generate random phases for 2x2 polarisation matrix elements
    phi = 360*rand(NumLinks,4,N_max,M);      % random phases for all users: [NumLinks pol path subpath]

    xpr_dB = randn(NumLinks,N_max,M);
    
    % get XPR distribution parameters (log-Normal)
    if LoSConnectionLinks
        xpr_mu_los = iterpar.LoS.xpr_mu;         % XPR mean [dB]
        xpr_sigma_los = iterpar.LoS.xpr_sigma;   % XPR std  [dB]
        xpr_dB(LoSConnectionLinks, 1:N(1),:)= xpr_dB(LoSConnectionLinks, 1:N(1),:)*xpr_sigma_los+xpr_mu_los;  % XPR [dB]
    end

    if NLoSConnectionLinks
        xpr_mu_nlos = iterpar.NLoS.xpr_mu;         % XPR mean [dB]
        xpr_sigma_nlos = iterpar.NLoS.xpr_sigma;   % XPR std  [dB]
        xpr_dB(NLoSConnectionLinks, 1:N(2),:)= xpr_dB(NLoSConnectionLinks, 1:N(2),:)*xpr_sigma_nlos+xpr_mu_nlos;  % XPR [dB]
    end

    % generate XPRs, dimensions are [NumLinks N M]
    xpr = 10.^(xpr_dB/10);     % XPR [linear]
    
    if N(1) < N(2)
        phi(LoSConnectionLinks,:,end+1-(N(2)-N(1)):end,:) = NaN;
        xpr(LoSConnectionLinks,end+1-(N(2)-N(1)):end,:) = NaN;
    elseif N(1) < N(2)
        phi(NLoSConnectionLinks,:,end+1-(N(1)-N(2)):end,:) = NaN;
        xpr(NLoSConnectionLinks,end+1-(N(1)-N(2)):end,:) = NaN;
    end

    % output
    bulk_parameters.subpath_phases = phi;        % in degrees
    bulk_parameters.xpr = xpr;                 % in linear scale

end
% end of output generation


% end of function stochastic





%%

%%%%%%%%%%%%%%%%%%% Variable definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to generate bulk parameters for B5 feeder scenarios
% See [1, Sec.5.3.1].
function bulk_parameters=static(wimpar,linkpar,antpar,fixpar,iterpar)

% NOTE! B5 scenario is only a CDL model with fixed parameters.

%First some checks.
if  ~(linkpar.MsVelocity==0) %KTH
    linkpar.MsVelocity=0;
end;%KTH

Scenario = iterpar.Scenario;
UserIndeces = iterpar.UserIndeces;
MsBsDistance = linkpar.MsBsDistance(UserIndeces);
NumLinks     = length(UserIndeces);

switch Scenario
    case {'B5a', 'B5b', 'B5c'}
        evalstr = sprintf('iterpar.LoS = fixpar.%s;',Scenario);
        eval(evalstr);
    case {'B5f'}
        evalstr = sprintf('iterpar.NLoS = fixpar.%s;',Scenario);
        eval(evalstr);
end

%extract line of sight probability
switch Scenario
    case {'B5a', 'B5b', 'B5c'}
        PropagCondition = ones(1,NumLinks);
    case {'B5f'}
        PropagCondition = zeros(1,NumLinks);
end

%indeces of LoS/NLoS links and the amount of them
iterpar.PropagCondition = PropagCondition;
LoSConnectionLinks = find(PropagCondition); iterpar.LoSConnectionLinks = LoSConnectionLinks;
NumLoSConnectionLinks = length(LoSConnectionLinks); iterpar.NumLoSConnectionLinks = NumLoSConnectionLinks;
NLoSConnectionLinks = find(PropagCondition==0); iterpar.NLoSConnectionLinks = NLoSConnectionLinks;
NumNLoSConnectionLinks = length(NLoSConnectionLinks); iterpar.NumNLoSConnectionLinks = NumNLoSConnectionLinks;

if strcmpi(Scenario,'b5b')%KTH
    if ~isfield(wimpar,'range')%KTH
        error('The field wimpar.range must exist in the stationary feeder scenario')%KTH
    end;%KTH
end;%KTH

if strcmpi(wimpar.AnsiC_core,'yes')%KTH
    error('The scenario B5 is only implemented in matlab not in ANSI-C.')%KTH
end;%KTH

%%%% Define local variables
M                     = wimpar.NumSubPathsPerPath;   %KTH
DelaySamplingInterval = wimpar.DelaySamplingInterval;%KTH

% make sure that user-specific parameters are row vectors
ThetaBs      = linkpar.ThetaBs(UserIndeces);
ThetaMs      = linkpar.ThetaMs(UserIndeces);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[taus_los,Pprime_los,Kcluster_los,...
    taus_nlos,Pprime_nlos,Kcluster_nlos] = fixedPdp(wimpar,iterpar);    % the same for each link

% Replace number of paths by number of tabulated paths
N(1) = length(find(~isnan(taus_los)));
N(2) = length(find(~isnan(taus_nlos)));
N_max = max(N(1),N(2));
iterpar.N = N;

taus_sorted = NaN*ones(NumLinks,N_max);
taus_sorted(LoSConnectionLinks,1:N(1)) = repmat(taus_los,NumLoSConnectionLinks,1);
taus_sorted(NLoSConnectionLinks,1:N(2)) = repmat(taus_nlos,NumNLoSConnectionLinks,1);

% Rounding to delay grid
if (DelaySamplingInterval>0)
    taus_rounded = DelaySamplingInterval*floor(taus_sorted/DelaySamplingInterval + 0.5);
else
    taus_rounded = taus_sorted;
end

P = NaN*ones(NumLinks,N_max);
P(LoSConnectionLinks,1:N(1)) = repmat(Pprime_los,NumLoSConnectionLinks,1);
P(NLoSConnectionLinks,1:N(2)) = repmat(Pprime_nlos,NumNLoSConnectionLinks,1);
P(LoSConnectionLinks,1:N(1)) = P(LoSConnectionLinks,1:N(1))./repmat(sum(P(LoSConnectionLinks,1:N(1)),2),1,N(1));
P(NLoSConnectionLinks,1:N(2)) = P(NLoSConnectionLinks,1:N(2))./repmat(sum(P(NLoSConnectionLinks,1:N(2)),2),1,N(2));


%transform Kcluster to linear
Kcluster = NaN*ones(NumLinks,1);
Kcluster(LoSConnectionLinks,:) = 10.^(Kcluster_los(1,:)/10);    % tranform to linear

offset = [0.0447 0.1413 0.2492 0.3715 0.5129 0.6797 0.8844 1.1481 1.5195 2.1551]; %KTH [1, Table 4-1] +/- offset angles, resulting Laplacian APS, with rms AS = 1 deg

% Determine AoDs %%
[AoD_path_los,iterpar.LoS.PerClusterAS_D,...
    AoD_path_nlos,iterpar.NLoS.PerClusterAS_D] = fixedAods(wimpar,iterpar);    % the same for each link

AoD_path = NaN*ones(NumLinks,N_max);
AoD_path(LoSConnectionLinks,1:N(1)) = repmat(AoD_path_los,NumLoSConnectionLinks,1);
AoD_path(NLoSConnectionLinks,1:N(2)) = repmat(AoD_path_nlos,NumNLoSConnectionLinks,1);

AoD_tmp = repmat(reshape(AoD_path.',1,N_max*NumLinks),M,1);

%Determine AoAs
[AoA_path_los,iterpar.LoS.PerClusterAS_A,...
    AoA_path_nlos,iterpar.NLoS.PerClusterAS_A] = fixedAoas(wimpar,iterpar);    % the same for each link

AoA_path = NaN*ones(NumLinks,N_max);
AoA_path(LoSConnectionLinks,1:N(1)) = repmat(AoA_path_los,NumLoSConnectionLinks,1);
AoA_path(NLoSConnectionLinks,1:N(2)) = repmat(AoA_path_nlos,NumNLoSConnectionLinks,1);

AoA_tmp = repmat(reshape(AoA_path.',1,N_max*NumLinks),M,1);

[offset_matrix_AoD, offset_matrix_AoA] = offset_matrix_generation(offset,iterpar);

%apply offset matrix
% NOTE! now array orientation parameter ThetaBs is disabled and AoD is always
% like in CDL model tables, 17.5.2006 PekKy.
theta_nm_aod    = AoD_tmp + offset_matrix_AoD;      % a (M x (NumLinks*N)) matrix
theta_nm_aoa    = AoA_tmp + offset_matrix_AoA;      % a (M x (NumLinks*N)) matrix

%Pair AoA rays randomly with AoD rays (within a cluster)
[dummy h]           = sort(rand(M,N_max*NumLinks),1);       % create N*NumLinks random permutations of integers [1:M]
inds                = h+repmat([1:M:M*N_max*NumLinks],M,1)-1;
theta_nm_aoa        = theta_nm_aoa(inds);    % rand

% Values of theta_nm_aoa and theta_nm_aod may be outside (-180,180).
% Wrapping of angles to range (-180,180)
theta_nm_aoa = prin_value(theta_nm_aoa);
theta_nm_aod = prin_value(theta_nm_aod);

%scatterer frequency
scatterer_freq = fixedScatterFreq(wimpar,iterpar);%KTH
scatterer_freq = repmat(scatterer_freq,1,NumLinks);%KTH

Phi_LOS = NaN*ones(NumLinks,1);
Phi_LOS(LoSConnectionLinks) = 360*(rand(NumLoSConnectionLinks,1)-0.5); %KTH
phi = 360*rand(NumLinks,N_max,M);        % random phases for all users

% employ the user-defined path loss model
[path_losses, linkpar, fixpar, iterpar] = feval(wimpar.PathLossModel,wimpar,linkpar,fixpar,iterpar);
path_losses = 10.^(-path_losses(:)/10);    % a (NumLinks x 1) vector
%path_losses = ones(NumLinks,1);

% Shadow-fading
% NOTE! all the links are fully uncorrelated, changed 12.12.2005 by Pekka
% BsNumber = linkpar.BsNumber(:);
% NumOfBs = max(BsNumber);
switch lower((Scenario))
    case {'b5a','b5c','b5f'}
        SF_sigma=3.4*ones(NumLinks,1);
    case 'b5b'
        lambda=3e8/wimpar.CenterFrequency;
        breakpoint_distance=4*(linkpar.MsHeight(UserIndeces)-1.6).*(linkpar.BsHeight(UserIndeces)-1.6)/lambda;
        within_breakpoint=MsBsDistance<breakpoint_distance;
        SF_sigma=(within_breakpoint*3+(~within_breakpoint)*7)';
end; %% switch
sigma_sf  = 10.^(0.1*(SF_sigma.*randn(NumLinks,1)));
%sigma_sf  = sigma_sf_all_sites(BsNumber);


% put AoDs, AoAs, scatter_freq and power?? gains into a 3D-array with dims [NumLinks N M]
theta_nm_aod=reshape(theta_nm_aod,M,N_max,NumLinks);  %KTH
theta_nm_aod=permute(theta_nm_aod,[3 2 1]);   %KTH
theta_nm_aoa=reshape(theta_nm_aoa,M,N_max,NumLinks);  %KTH
theta_nm_aoa=permute(theta_nm_aoa,[3 2 1]);   %KTH
scatterer_freq=reshape(scatterer_freq,M,N_max,NumLinks);  %KTH
scatterer_freq=permute(scatterer_freq,[3 2 1]);   %KTH

bulk_parameters=struct( 'delays',taus_rounded,...
    'path_powers',P,...                 %before: 'subpath_powers',Psub,...
    'aods',theta_nm_aod,...             %
    'aoas',theta_nm_aoa,...
    'subpath_phases',phi,...
    'Kcluster',Kcluster,...           % in dB.
    'Phi_LOS',Phi_LOS,...               % phases for LOS paths, in degrees
    'path_losses',path_losses,...       % in linear scale
    'shadow_fading',sigma_sf,...          % in linear scale
    'MsBsDistance',MsBsDistance,... % This output is needed since the originally generated MsBsDistances are fitted inside the applicabity ranges of the Scenarios
    'scatterer_freq',scatterer_freq,...
    'propag_condition', PropagCondition.',...
    'LoSO2ILinks', linkpar.LoSO2ILinks,...
    'LoSO2VLinks', linkpar.LoSO2VLinks,...
    'NLoSO2ILinks', linkpar.NLoSO2ILinks,...
    'NLoSO2VLinks', linkpar.NLoSO2VLinks,...
    'o2v_shadow_fading', linkpar.sigma_sf_o2v,...
    'user_indeces',UserIndeces.');%KTH

%-- STEP 9 - STEP10b --%
if strcmpi(wimpar.PolarisedArrays,'yes')

    % generate random phases for 2x2 polarisation matrix elements
    phi = 360*rand(NumLinks,4,N_max,M);      % random phases for all users: [NumLinks pol path subpath]

    xpr_dB = randn(NumLinks,N_max,M);

    % get XPR distribution parameters (log-Normal)
    if LoSConnectionLinks
        xpr_mu_los = iterpar.LoS.xpr_mu;         % XPR mean [dB]
        xpr_sigma_los = iterpar.LoS.xpr_sigma;   % XPR std  [dB]
        xpr_dB(LoSConnectionLinks, 1:N(1),:)= xpr_dB(LoSConnectionLinks, 1:N(1),:)*xpr_sigma_los+xpr_mu_los;  % XPR [dB]
    end

    if NLoSConnectionLinks
        xpr_mu_nlos = iterpar.NLoS.xpr_mu;         % XPR mean [dB]
        xpr_sigma_nlos = iterpar.NLoS.xpr_sigma;   % XPR std  [dB]
        xpr_dB(NLoSConnectionLinks, 1:N(2),:)= xpr_dB(NLoSConnectionLinks, 1:N(2),:)*xpr_sigma_nlos+xpr_mu_nlos;  % XPR [dB]
    end

    % generate XPRs, dimensions are [NumLinks N M]
    xpr = 10.^(xpr_dB/10);     % XPR [linear]

    if N(1) < N(2)
        phi(LoSConnectionLinks,:,end+1-(N(2)-N(1)):end,:) = NaN;
        xpr(LoSConnectionLinks,end+1-(N(2)-N(1)):end,:) = NaN;
    elseif N(1) < N(2)
        phi(NLoSConnectionLinks,:,end+1-(N(1)-N(2)):end,:) = NaN;
        xpr(NLoSConnectionLinks,end+1-(N(1)-N(2)):end,:) = NaN;
    end

    % output
    bulk_parameters.subpath_phases = phi;        % in degrees
    bulk_parameters.xpr = xpr;                 % in linear scale

end
%%

