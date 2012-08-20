%WIM_CORE Channel coefficient computation for a geometric channel model
%   [H DELTA_T FINAL_PHASES FINAL_PHASES_LOS]=WIM_CORE(WIMPAR,LINKPAR,ANTPAR,BULKPAR,BSGAIN,BSGAIN_LOS,MSGAIN,MSGAIN_LOS,OFFSET_TIME, BSGAINISSCALAR, MSGAINISSCALAR, PCind)
%   This is the wim_core aka the big for loop. It implements the formula in
%   [2, Eq. 3.26-28].
%
%   Outputs:
%
%   H               - [UxSxNxTxK] array of channel coefficients
%   DELTA_T         - time sampling intervals (in seconds) for all links
%   FINAL_PHASES    - final phases of all subpaths in degrees over (-180,180)
%   FINAL_PHASES_LOS- final phases for LOS paths in degrees over (-180,180)
%   
%   Inputs:
%
%   WIMPAR          - input struct, see WIMPARSET
%   LINKPAR         - input struct, see LINKPARSET
%   ANTPAR          - input struct, see ANTPARSET
%   BULKPAR         - input BULKPAR, see GENERATE_BULK_PAR
%   BSGAIN          - [KxSxNxM] array of interpolated antenna field
%                     patterns (complex)
%   BSGAIN_LOS      - [KxS] array of interpolated antenna field patterns
%                     (complex) for LOS paths. Only used with the LOS
%                     option; it is set to scalar otherwise.
%   MSGAIN          - [KxUxNxM] array of interpolated antenna field
%                     patterns (complex)
%   MSGAIN_LOS      - [KxU] array of interpolated antenna field patterns
%                     (complex) for LOS paths. Only used with the LOS
%                     option; it is set to scalar otherwise.
%   OFFSET_TIME     - time offset added to the initial phase (set to zero by default)
%   BSGAINISSCALAR  - this is 1 if BsGain is uniform over azimuth, 0 otherwise.
%   MSGAINISSCALAR  - this is 1 if MsGain is uniform over azimuth, 0 otherwise.
%
%   With 'polarized' option:
%
%   BSGAIN          - [KxSx2xNxM] array of interpolated antenna field
%                     patterns (complex), where the third dimension are
%                     the patterns for [V H] polarizations. 
%   MSGAIN          - [KxUx2xNxM] array of interpolated antenna field
%                     patterns (complex), where the third dimension are
%                     the patterns for [V H] polarizations. 
%
%   The ANCI-C version is develeped in WINNER project originally to
%   implement the 3GPP SCM.
%   To compile the ANSI-C written optimized core, type
%
%       mex scm_mex_core.c
%
%   at MATLAB prompt. For further documentation on the ANSI-C implementation 
%   of the WIM_CORE, see SCM_MEX_CORE.
%
%   Ref. [1]: 3GPP TR 25.996 v6.1.0 (2003-09)
%        [2]: D1.1.2 V1.0, "WINNER II channel models"
%
%   Authors: Giovanni Del Galdo (TUI), Marko Milojevic (TUI), Jussi Salmi (HUT),  
%   Christian Schneider (TUI), Jari Salo (HUT), Pekka Ky�sti (EBIT), 
%   Daniela Laselva (EBIT), Lassi Hentil� (EBIT)

% Bug fixes:
%   Erroneus variable name s changed to u on line 573. Caused error with
%    asymmetric MIMO configuration (e.g.2x4) & CDL models or B5 scenario. (22.8.2006 PekKy)
%   Polarised arrays, non-ANSI-C version, power scaling changed harmonised
%    with ANSI-C version on rows 419-420, 462-468. (22.8.2006 PekKy)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   --------                     %%
function [H, delta_t, output_SubPathPhases, output_Phi_LOS] = wim_core (wimpar,linkpar,antpar,bulkpar,BsGain,BsGain_Theta_BS,MsGain,MsGain_Theta_MS,offset_time, BsGainIsScalar, MsGainIsScalar,PCind)
%%                   --------                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% offset_time [samples] = defines the starting point (in samples) of the
%                         time axis
% Examples: you want to calculate 1000 time samples calling the wim_core
%           twice (everytime for 500 timesamples)
%           H1 = wim_core (..., 0)
%           H2 = wim_core (..., 500)
%



DEBUG_MODE_FLAG   = 0;
PROFILE_MODE_FLAG = 0;
DISPLAY_MODE_FLAG = 0;

S=size(BsGain,2);                     % number of receiving antennas
U    = size(MsGain,2);                % number of transmitting antennas
N    = size(bulkpar.delays,2);        % number of paths
T    = wimpar.NumTimeSamples;         % number of time samples
K    = length(PCind);                 % number of links
M    = wimpar.NumSubPathsPerPath;     % number of subpaths

AnsiC_core = wimpar.AnsiC_core;


% intra-cluster delays spread added based on the D1.1.1
if strcmpi(wimpar.IntraClusterDsUsed,'yes')
    AnsiC_core = 'no';      % NOTE! IntraClusterDsUsed forces AnsiC_core = 'no', change this after scm_mex_core.c is modified
    
    NumRaysPerSubCluster = [10,6,4];
    RayOrder = [1,2,3,4,5,6,7,8,19,20,9,10,11,12,17,18,13,14,15,16];
    P = bulkpar.path_powers; P(isnan(P))=-Inf;
    SortedPower = fliplr(sort(P,2)); P(isinf(P))=NaN;
    for xx = 1:size(P,1)
        SubClusterInd(xx,:) = P(xx,:) > SortedPower(xx,3); % Index of the cluster to be divided
    end
    
else    % special case, one midpath only
    LM = 1:M;
    LN=M;
    L = length(LN);
    P = bulkpar.path_powers;
    SubClusterInd = zeros(size(P));
end



H = zeros(U,S,N+4,T,K); %"+4" is due to 4 extra sub-clusters


% define element spacing vectors if scalars are given
if (length(antpar.MsElementPosition)==1)
    antpar.MsElementPosition=[0:antpar.MsElementPosition:antpar.MsElementPosition*(U-1)];
end

if (length(antpar.BsElementPosition)==1)
    antpar.BsElementPosition=[0:antpar.BsElementPosition:antpar.BsElementPosition*(S-1)];
end


% Set internal parameters
speed_of_light=2.99792458e8;
wavelength=speed_of_light/wimpar.CenterFrequency;

% dummy
output_Phi_LOS       = zeros(K,1);

% let's make the time axis - for that we need to check UniformTimeSampling
% and the MSs' velocities
% Note: SampleDensity is samples per half wavelength.
if strcmp(wimpar.UniformTimeSampling,'yes')
    
    max_vel = max(linkpar.MsVelocity);
    delta_t = repmat((wavelength / max_vel)/2/wimpar.SampleDensity,K,1);
    
else % 'UniformTimeSampling' is 'no'
    
    delta_t = (wavelength ./ linkpar.MsVelocity(PCind).')./2/wimpar.SampleDensity ;
    
end
t = repmat(delta_t,1,T).*repmat([0:T-1]+offset_time,K,1); % matrix containing the time axes for all links [KxT]

% Time axis generation for fixed feeder links (B5 scenarios)
tmp = zeros(1,length(linkpar.ScenarioVector));
tmp(PCind) = 1;
B5ind = find((linkpar.ScenarioVector>=7 & linkpar.ScenarioVector<=9).*tmp);

if length(B5ind)>0
    H = zeros(U,S,N+4,T,K);  %"+4" is due to 4 extra sub-clusters (even though not used with B5)
    SubClusterInd(B5ind,:) = zeros(length(B5ind),size(bulkpar.delays,2));
    AnsiC_core = 'no';      % NOTE! B5 forces AnsiC_core = 'no', change this after scm_mex_core.c is modified
    for k=1:length(B5ind) 
        B5ind2(k) = find((B5ind(k)==PCind));
    end
    % not final
    wimpar.TimeVector=linspace(0,wimpar.end_time,T);
    % not final
    
    t(B5ind2,:) = repmat(wimpar.TimeVector,length(B5ind),1);%KTH
    linkpar.MsVelocity(B5ind)=zeros(length(B5ind),1);
    delta_t(B5ind2)=repmat(wimpar.TimeVector(2)-wimpar.TimeVector(1),length(B5ind),1); %% Dummy value
end

k_CONST = 2*pi/wavelength;      % wave number


%%%%%%%%%%%%%%%%%%ANSI-C core part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check if ANSI-C core is used
if strcmpi(wimpar.AnsiC_core,'yes')
    
    % different modes
    GENERAL = 1;
    POLARIZED = 2;
    LOS = 3;   
    
    look_up_points = wimpar.LookUpTable; % set if look-up table is used for sin/cos
        
    if (BsGainIsScalar && MsGainIsScalar)
        GainsAreScalar = 1;
    else 
        GainsAreScalar = 0;
    end
    
    
    if ~strcmpi(wimpar.PolarisedArrays,'yes')
        
        if DISPLAY_MODE_FLAG
            disp('entering main loop...');
        end
        
        % adjusting parameters for calling the C routine
        d_u = antpar.MsElementPosition*wavelength;
        d_s = antpar.BsElementPosition*wavelength;
        aod = bulkpar.aods(PCind,1:N,1:M)*pi/180;
        aoa = bulkpar.aoas(PCind,1:N,1:M)*pi/180;  
        phase = bulkpar.subpath_phases(PCind,1:N,1:M)*pi/180;
        v = linkpar.MsVelocity(PCind);
        theta_v = linkpar.MsDirection(PCind)*pi/180;
        sq_Pn = sqrt(bulkpar.path_powers(PCind,1:N*L));
        
        % calling the C-mex routine for general coefficients 
        [H output_SubPathPhases] = scm_mex_core(GENERAL, BsGain(PCind,:,:,:), MsGain(PCind,:,:,:), aod, aoa,...
            d_s, d_u, phase, t, k_CONST, v, theta_v, sq_Pn, look_up_points, U, S, N, L, M, K, T, GainsAreScalar, LM, LN);
        
        %output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases)); %changed due tests
        output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases(1:K,1:N,1:M)));
        
    else % it's polarized!
        if DISPLAY_MODE_FLAG
            disp('entering polarized option...');
        end    
        
        output_SubPathPhases = zeros(K,4,N,M);
        %temp_output_SubPathPhases = zeros(K,N,M);
        
        % BsGain must have size: [K x S x 2 x N x M]
        % the first dimension in the polarization must be vertical
        %
        % MsGain must have size: [K x U x 2 x N x M]
        % the first dimension in the polarization must be vertical
        %
        % subpath_phases has size: [K x 4 x N x M]
        %
        % bulkpar.xpd has size: [K,2,N]
        % adjusting parameters for calling the C routine
        d_u = antpar.MsElementPosition * wavelength;
        d_s = antpar.BsElementPosition * wavelength;
        
        X_BS_v = reshape(BsGain(:,:,1,:,:),K,S,N,M); 
        X_BS_h = reshape(BsGain(:,:,2,:,:),K,S,N,M); 
        X_MS_v = reshape(MsGain(:,:,1,:,:),K,U,N,M); 
        X_MS_h = reshape(MsGain(:,:,2,:,:),K,U,N,M);     
        aod = bulkpar.aods(PCind,1:N,1:M)*pi/180;
        aoa = bulkpar.aoas(PCind,1:N,1:M)*pi/180;  
        phase_v_v = reshape(bulkpar.subpath_phases(PCind,1,1:N,1:M),K,N,M)*pi/180;
        phase_v_h = reshape(bulkpar.subpath_phases(PCind,2,1:N,1:M),K,N,M)*pi/180;
        phase_h_v = reshape(bulkpar.subpath_phases(PCind,3,1:N,1:M),K,N,M)*pi/180;
        phase_h_h = reshape(bulkpar.subpath_phases(PCind,4,1:N,1:M),K,N,M)*pi/180;
        %sq_r_n1 = reshape(sqrt(1/bulkpar.xpd(1:K,1,1:N)),K,N);
        %sq_r_n2 = reshape(sqrt(1/bulkpar.xpd(1:K,2,1:N)),K,N);        
%         r_n1 = 1 ./ reshape(bulkpar.xpd(1:K,1,1:N),K,N); % Jari Apr 17, 2005
%         r_n2 = 1 ./ reshape(bulkpar.xpd(1:K,2,1:N),K,N); % Jari Apr 17, 2005   
        r_n1 = 1 ./ bulkpar.xprV(PCind,1:N,1:M); 
        r_n2 = 1 ./ bulkpar.xprH(PCind,1:N,1:M);   

        v = linkpar.MsVelocity(PCind);
        theta_v = linkpar.MsDirection(PCind)*pi/180;
        sq_Pn = sqrt(bulkpar.path_powers(PCind,1:N*L));
        
        % the ANSI-C function call
        [H temp_output_SubPathPhases] = scm_mex_core(POLARIZED, X_BS_v, X_BS_h, X_MS_v, X_MS_h, aod, aoa, d_s, d_u, phase_v_v, phase_v_h, phase_h_v, phase_h_h, r_n1, r_n2, t, k_CONST, v, theta_v, sq_Pn, look_up_points, U, S, N, L, M, K, T, GainsAreScalar, LM, LN);
        
        output_SubPathPhases(:,1,:,:) = temp_output_SubPathPhases;
        output_SubPathPhases(:,2,:,:) = temp_output_SubPathPhases;
        output_SubPathPhases(:,3,:,:) = temp_output_SubPathPhases;
        output_SubPathPhases(:,4,:,:) = temp_output_SubPathPhases;
        
        output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases));
        
        
        
    end % is it polarized?
    
    %%%%%%
    if PROFILE_MODE_FLAG
        profile report
    end
    %%%%%%
    
    % LOS OPTION
    
    if (bulkpar.propag_condition(PCind(1))==1)      % if LOS links
        
        
        % Take the values of K factors and probability of having LOS case
        K_factors       = bulkpar.Kcluster;%;K_factors;
        
        %indx            = (K_factors.'~=0)';
        
        ThetaBs      = linkpar.ThetaBs; ThetaBs=ThetaBs(:).';
        ThetaMs      = linkpar.ThetaMs; ThetaMs=ThetaMs(:).';
        
        % if strcmp(str,'safe') % we only do it 'safe' for this option
        
        %adjusting parameters for calling the C-language routine
        output_Phi_LOS       = zeros(length(linkpar.ScenarioVector),1);
        d_u = antpar.MsElementPosition * wavelength;
        d_s = antpar.BsElementPosition * wavelength;
        
        
        % the ANSI-C function call
        [H output_Phi_LOS] = scm_mex_core(LOS, BsGain_Theta_BS(PCind,:), MsGain_Theta_MS(PCind,:), ThetaBs*pi/180, ThetaMs*pi/180, d_s, d_u, bulkpar.Phi_LOS(:,1)* pi/180, t, k_CONST, linkpar.MsVelocity(:), linkpar.MsDirection(1,:)*pi/180, H, output_Phi_LOS, K_factors, U, S, N*L, K, T);
        
        % adjusting angles
        output_Phi_LOS = prin_value((output_Phi_LOS*180/pi + bulkpar.Phi_LOS));
        
    end % if 'LOS' condition   
    %%%%%%%%%%%%%%%%%%ANSI-C core part ends%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
else % if ANSI-C is not used
    if ~strcmpi(wimpar.PolarisedArrays,'yes')
        
        if DISPLAY_MODE_FLAG
            disp('entering main loop...');
        end
        
        output_SubPathPhases = zeros(K,N,M);
        
        for u = 1:U % cycles (MS) antennas
            %du = antpar.MsElementSpacingULA * (u-1) * wavelength;
            du = antpar.MsElementPosition(u)*wavelength;
            for s = 1:S % cycles Tx (BS) atennas
                %ds = antpar.BsElementSpacingULA * (s-1) * wavelength;
                ds = antpar.BsElementPosition(s)*wavelength;
                for kk = 1:K % cycles links
                    k = PCind(kk);
                    LN_index = 0; %
                    for n = 1:N % cycles paths
                        
                        LM_index = 0; %
                         
                        
                        if SubClusterInd(kk,n)==1; L=3; LM = RayOrder; LN = NumRaysPerSubCluster;
                           path_powers = [P(kk,n)*10/20 P(kk,n)*6/20 P(kk,n)*4/20]; 
                        else L=1; LN=M; LM=1:M;
                            path_powers = P(kk,n);
                        end
                        
                        for km = 1:L % cycles subclusters 
                            
                            LN_index = LN_index+1; %Running index of the clusters (including sub-clusters)
                            
                            temp = zeros(M,T);   % oversized, just to keep it always the same size
                            
                            for m=1:LN(km) % cycles rays
                                
                                LM_index = LM_index+1; % Running index of the rays
                                
                                % Calculate Doppler frequency nu of scatterer m
                                if sum(k==B5ind)    % IF current link is B5
                                    nu = bulkpar.scatterer_freq(k,n,m);
                                else    % IF not B5
                                    nu = (linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,LM(LM_index)) - linkpar.MsDirection(k))*pi/180))/wavelength;
                                end
                                
                                temp(m,:) =  BsGain(k,s,n,LM(LM_index)) *...
                                    exp(j*(...
                                    k_CONST * ds * sin((bulkpar.aods(k,n,LM(LM_index)))*pi/180) +...
                                    (bulkpar.subpath_phases(k,n,LM(LM_index))*pi/180)+...
                                    k_CONST * du * sin((bulkpar.aoas(k,n,LM(LM_index)))*pi/180)...
                                    )) *...
                                    MsGain(k,u,n,LM(LM_index)) * exp(1j*2*pi*nu * t(kk,:) );
                                
                                
                            end % rays

                            H(u,s,LN_index,:,kk) =  sqrt(path_powers(km) / LN(km)) * sum(temp,1);
                            
                        end % subclusters
                        
                    end % paths 
                    
                end % links 
            end % Tx antennas 
        end % Rx antennas 
        
        
        for kk = 1:K % cycles links  % Of course the for loop could be avoided. 
            k = PCind(kk);
            for n = 1:N % cycles paths
                for m=1:M % cycles supaths
                    
                    % SIMPLE LOOP
                    output_SubPathPhases(kk,n,m) =  k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,m) - linkpar.MsDirection(k))*pi/180) * (delta_t(kk)*T);   
                    
                end % subpaths
            end % paths 
        end % links 
        
        output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases(PCind,:,:)));
        
        
        
    else % it's polarized!
        if DISPLAY_MODE_FLAG
            disp('entering polarized option...');
        end    
        
        output_SubPathPhases = zeros(K,4,N,M);
        
        % Set polarisation matrix powers according to XPRs
        % Assume Pvv+Pvh=Phh+Phv=1. In this case
%         Pvv = bulkpar.xprV./(1+bulkpar.xprV);
%         Pvh = 1-bulkpar.xprV./(1+bulkpar.xprV);
%         Phh = bulkpar.xprH./(1+bulkpar.xprH);
%         Phv = 1-bulkpar.xprH./(1+bulkpar.xprH);
        % Pxy replaced by R_ni (8.5.2006 PekKy)
        r_n1 = 1 ./ bulkpar.xpr;
        r_n2 = r_n1; % D1.1.2 definition
        %r_n2 = 1 ./ bulkpar.xprH; 
        
        
        % BsGain has must have size: [K x S x 2 x N x M]
        % the first dimension in the polarization must be vertical
        %
        % MsGain has must have size: [K x U x 2 x N x M]
        % the first dimension in the polarization must be vertical
        %
        % subpath_phases has size: [K x 4 x N x M]
        % bulkpar.xpd has size: [K,2,N]
        temp = zeros(M,T);    
        for u = 1:U % cycles (MS) antennas
            %du = antpar.MsElementSpacingULA * (u-1) * wavelength;
            du = antpar.MsElementPosition(u)*wavelength;
            for s = 1:S % cycles Tx (BS) atennas
                %ds = antpar.BsElementSpacingULA * (s-1) * wavelength;
                ds = antpar.BsElementPosition(s)*wavelength;
                for kk = 1:K % cycles links
                    k = PCind(kk);
                    LN_index = 0; %
                    for n = 1:N % cycles paths
                        
                        LM_index = 0; % 
                        
                        if SubClusterInd(kk,n)==1; L=3; LM = RayOrder; LN = NumRaysPerSubCluster;
                           path_powers = [P(kk,n)*10/20 P(kk,n)*6/20 P(kk,n)*4/20]; 
                        else L=1; LN=M; LM=1:M;
                            path_powers = P(kk,n);
                        end
                        
                        for km = 1:L % cycles midpaths
                            
                            LN_index = LN_index+1; %Running index of the clusters (including sub-clusters)
                            
                            temp = zeros(M,T);   % oversized, just to keep it always the same size
                            
                            for m=1:LN(km) % cycles subpaths
                                
                                LM_index = LM_index+1;
                                                              
%                                 % Assume Pvv+Pvh=Phh+Phv=1. In this case        % Commented 22.8.2006, PekKy  
%                                 % Pvv=xprV/(1+xprV) and Pvh=1-xprV/(1+xprV),
%                                 % Phh=xprH/(1+xprH) and Phv=1-xprH/(1+xprH)
%                                 temp(m,:) =  [BsGain(k,s,1,n,LM(LM_index)) BsGain(k,s,2,n,LM(LM_index))] * ...
%                                     [ sqrt(Pvv(k,n,m)) * exp(j*bulkpar.subpath_phases(k,1,n,LM(LM_index))*pi/180),   sqrt(Pvh(k,n,m)) * exp(j*bulkpar.subpath_phases(k,2,n,LM(LM_index))*pi/180);...
%                                       sqrt(Phv(k,n,m)) * exp(j*bulkpar.subpath_phases(k,3,n,LM(LM_index))*pi/180),   sqrt(Phh(k,n,m)) * exp(j*bulkpar.subpath_phases(k,4,n,LM(LM_index))*pi/180)] *...
%                                     [MsGain(k,u,1,n,LM(LM_index)); MsGain(k,u,2,n,LM(LM_index))] * ...
%                                     exp(j * (k_CONST * ds * sin((bulkpar.aods(k,n,LM(LM_index)))*pi/180))) * ...
%                                     exp(j * (k_CONST * du * sin((bulkpar.aoas(k,n,LM(LM_index)))*pi/180))) * ...                    
%                                     exp(j * k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,LM(LM_index)) - linkpar.MsDirection(k))*pi/180) * t(k,:));

                                % Calculate Doppler frequency nu of scatterer m
                                if sum(k==B5ind)    % IF current link is B5
                                    nu = bulkpar.scatterer_freq(k,n,m);
                                else    % IF not B5
                                    nu = (linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,LM(LM_index)) - linkpar.MsDirection(k))*pi/180))/wavelength;
                                end

                                temp(m,:) =  [BsGain(k,s,1,n,LM(LM_index)) BsGain(k,s,2,n,LM(LM_index))] * ...
                                    [ exp(j*bulkpar.subpath_phases(k,1,n,LM(LM_index))*pi/180)  ,  sqrt(r_n1(k,n,m)) * exp(j*bulkpar.subpath_phases(k,2,n,LM(LM_index))*pi/180);...
                                    sqrt(r_n2(k,n,m)) * exp(j*bulkpar.subpath_phases(k,3,n,LM(LM_index))*pi/180)  ,  exp(j*bulkpar.subpath_phases(k,4,n,LM(LM_index))*pi/180)] *...
                                    [MsGain(k,u,1,n,LM(LM_index)); MsGain(k,u,2,n,LM(LM_index))] * ...
                                    exp(j * (k_CONST * ds * sin((bulkpar.aods(k,n,LM(LM_index)))*pi/180))) * ...
                                    exp(j * (k_CONST * du * sin((bulkpar.aoas(k,n,LM(LM_index)))*pi/180))) * ...   
                                    exp(1j*2*pi*nu * t(kk,:) );
                                
                                
                            end % rays
                            
                            H(u,s,LN_index,:,kk) =  sqrt(path_powers(km) / LN(km)) * sum(temp,1);
                            
                        end % subclusters
                    end % paths 
                end % links 
            end % Tx antennas 
        end % Rx antennas 
        
%         for k = 1:K % cycles links  % Of course the for loop could be avoided. 
%             for n = 1:N % cycles paths
%                 for m=1:M % cycles supaths
%                     
%                     
%                     output_SubPathPhases(k,:,n,m) =  (k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,m) - linkpar.MsDirection(k))*pi/180) * (t(k,end)+delta_t(k))) * ones(1,4);
%                     
%                 end % subpaths
%             end % paths 
%         end % links 
%         
%         output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases));
        
        
         for kk = 1:K % cycles links  % Of course the for loop could be avoided. 
            k = PCind(kk);
            for n = 1:N % cycles paths
                for m=1:M % cycles supaths
                    
                    % SIMPLE LOOP
                    output_SubPathPhases(kk,:,n,m) =  (k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,m) - linkpar.MsDirection(k))*pi/180) * (delta_t(kk)*T)) * ones(1,4);
                    
                end % subpaths
            end % paths 
        end % links 
        
        output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases(PCind,:,:,:)));
        
        
        
    end % is it polarized?
    
    
    
    
    %%% LOS OPTION %%%
    
    % index to LOS but not B5 links
    LosNonB5ind = find(bulkpar.propag_condition'.*(linkpar.ScenarioVector<7 | linkpar.ScenarioVector>9));
    
    % If LOS, but not B5 link and not 'FixedPDP'
    if (bulkpar.propag_condition(PCind(1))==1) & length(LosNonB5ind)>0    
        
        % Take the values of K factors
        K_factors       = bulkpar.Kcluster;%;K_factors;
        
        ThetaBs      = linkpar.ThetaBs; ThetaBs=ThetaBs(:).';
        ThetaMs      = linkpar.ThetaMs; ThetaMs=ThetaMs(:).';
        
        output_Phi_LOS       = zeros(length(linkpar.ScenarioVector),1);
        
        for kk = 1:length(LosNonB5ind) % cycles links
            k = LosNonB5ind(kk);
            k_ind = find(k==PCind);     % index to current LOS/NLOS links which are LOS but not B5
            output_Phi_LOS(k,1) = k_CONST * linkpar.MsVelocity(k) * cos((ThetaMs(k) - linkpar.MsDirection(k))*pi/180) * (t(k_ind,end)+delta_t(k_ind));
            for u = 1:U % cycles (MS) antennas
                du = antpar.MsElementPosition(u)*wavelength;
                for s = 1:S % cycles (BS) antennas
                    ds = antpar.BsElementPosition(s)*wavelength;
                    temp =  BsGain_Theta_BS(k,s) * exp(j * k_CONST * ds * sin( ThetaBs(k)*pi/180)).* ...
                            MsGain_Theta_MS(k,u) * exp(j * (k_CONST * du * sin( ThetaMs(k)*pi/180  ) + bulkpar.Phi_LOS(k,1) * pi/180 )) * ...
                            exp(j * k_CONST * linkpar.MsVelocity(k) * cos((ThetaMs(k) - linkpar.MsDirection(k))*pi/180) * t(find(k==PCind),:));
                        
                    
                    H(u,s,1,:,k_ind)= (sqrt(1/(K_factors(k)+1)) * squeeze(H(u,s,1,:,k_ind)) + sqrt(K_factors(k)/(K_factors(k)+1)) * temp.').';
                end % Rx antennas
                
            end % Tx antennas
            
            H(:,:,2:end,:,k_ind)= sqrt(1/(K_factors(k)+1)) *  H(:,:,2:end,:,k_ind);      
        end % links
        
        output_Phi_LOS = prin_value((output_Phi_LOS*180/pi + bulkpar.Phi_LOS));
    end % if 'LOS' propagation condition
    
end  % end if ANSI-C is used


%% B5 links %%%
%   direct rays are added and cluster powers adjusted
%   according to cluster-wise K-factors given in [2, tables 7.17-24]
%
% index to links which are LOS and B5
LosB5ind = find(bulkpar.propag_condition'.*(linkpar.ScenarioVector>=7 & linkpar.ScenarioVector<=9));

if (bulkpar.propag_condition(PCind(1))==1) & length(LosB5ind)>0 
    
    Kcluster = bulkpar.Kcluster;             % read cluster-wise K-factors
    
%    output_Phi_LOS       = zeros(K,1);
    
    for kk = 1:length(LosB5ind) % cycles links
        k = LosB5ind(kk);
        k_ind = find(k==PCind);     % index to current LOS/NLOS links which are LOS and B5
        for u = 1:U % cycles (MS) antennas
            du = antpar.MsElementPosition(u)*wavelength;
            for s = 1:S % cycles (BS) atennas
                
                n = 1;     % index to cluster with a direct ray, in WIM2 always 1st cluster
                ds = antpar.BsElementPosition(s)*wavelength;
                
                aod_direct = bulkpar.aods(k,n,1)-bulkpar.aods(k,n,2);     % AoD for the direct ray (middle)
                aoa_direct = bulkpar.aoas(k,n,1)-bulkpar.aoas(k,n,2);     % AoA for the direct ray (middle)
                
                % antenna gain of direct ray is approximated by linear interpolation 
                BsGain_direct = mean(BsGain(k,s,n,1:2));    
                MsGain_direct = mean(MsGain(k,u,n,1:2));        % 22.8.2006 PekKy, index s corrected to u
                
                nu = 0;     % LOS ray has always 0 Hz Doppler in B5 scenarios
                
                temp =  BsGain_direct * exp(j * k_CONST * ds * sin( aod_direct*pi/180))* ...
                        MsGain_direct * exp(j * (k_CONST * du * sin( aoa_direct*pi/180  ) + bulkpar.Phi_LOS(k,1) * pi/180 )) * ...
                        exp(1j*2*pi*nu * t(k_ind,:) );
                
                H(u,s,n,:,k_ind) = (sqrt(1/(Kcluster(k,1)+1)) * squeeze(H(u,s,n,:,k_ind)) + sqrt(Kcluster(k,1)/(Kcluster(k,1)+1)) * temp.').';
                
                output_Phi_LOS(k_ind,n) = 0;
                
            end % Rx antennas
        end % Tx antennas
             
    end % links
    
    output_Phi_LOS(k_ind,:) = prin_value((output_Phi_LOS(k_ind,:)*180/pi + bulkpar.Phi_LOS(k_ind,:)));

end     % end of B5 part

%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%
%%%%%        That's all folks !!!
%%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);