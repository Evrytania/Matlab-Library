function [H, delays, full_output]=wim(wimpar,linkpar,antpar,initvalues)
%WIM  WINNER Phase II Channel Model (WIM2)
%   H=WIM(WIMPAR,LINKPAR,ANTPAR) is a 5D-array of channel coefficients. For
%   explanation of the input parameter structs, see WIMPARSET, LINKPARSET,
%   and ANTPARSET. SIZE(H)=[U S N T K], where U is the number of MS (RX)
%   elements, S is the number of BS (TX) elements, N is the number of paths,
%   T is the number of time samples, and K is the number of links. If K=1,
%   the final dimension will be dropped, i.e. H is a 4D-array.
%
%   [H DELAYS]=WIM(...) outputs also a [KxN] matrix of path delays (in seconds).
%
%   [H DELAYS BULKPAR]=WIM(...) outputs also the struct BULKPAR, whose fields
%   are as follows:
%
%   With NLOS propagation condition:
%
%   delays          - path delays in seconds [KxN]
%   path_powers     - relative path powers [KxN]
%   aods            - angles of departure in degrees over (-180,180) [KxNxM]
%   aoas            - angles of arrival in degrees over (-180,180) [KxNxM]
%   subpath_phases  - final phases for subpaths in degrees over (0,360) [KxNxM]
%   path_losses     - path losses in linear scale [Kx1]
%   MsBsDistance    - distances between MSs and BSs in meters [1xK]
%   shadow_fading   - shadow fading losses in linear scale [Kx1]
%   delta_t         - time sampling intervals for all links [Kx1]
%
%   In addition, when LOS condition (in addition to the above):
%
%   K_factors       - K factors for all links [Kx1]
%   Phi_LOS         - final phases for LOS paths in degrees over (-180,180) [Kx1]
%
%   [H ...]=WIM(...,INIT_VALUES) uses initial values given in the struct
%   INIT_VALUES, instead of random parameter generation. INIT_VALUES has
%   the same format as BULKPAR, except that SUBPATH_PHASES are now the
%   initial phases. Also, time sampling intervals (delta_t) are not used
%   (they are recalculated for every call of WIM).
%
%   Example:
%       % to generate matrices for 10 links with default parameters
%       H=wim(wimparset,linkparset(10),antparset);
%       % to generate matrices for A1 LOS scenario
%       wimpar=wimparset;wimpar.Scenario='A1';wimpar.PropagCondition: 'los'
%       H=wim(wimpar,linkparset(10),antparset);
%
%   Ref. [1]: D1.1.2 V1.0, "WINNER II channel models"
%        [2]: 3GPP TR 25.996 v6.1.0 (2003-09)
%
%   See also WIMPARSET, LINKPARSET, ANTPARSET

%   Authors: Jari Salo (HUT), Giovanni Del Galdo (TUI), Pekka Ky�sti (EBIT),
%   Daniela Laselva (EBIT), Marko Milojevic (TUI), Christian Schneider (TUI)
%   Lassi Hentil� (EBIT), Mikko Alatossava (CWC/UOULU),Zhiwen Wu(BUPT),Yu
%   Zhang(BUPT),Jianhua Zhang(BUPT),Guangyi Liu(CMCC)


% Note: all units are in degrees, meters, Hertz (1/s) and meters/second (m/s)




ni=nargin;
if (ni<3 || ni>4)
    error('WIM requires three or four input arguments !')
end


% Read fixed scenario dependent parameters from a table
fixpar = ScenParTables(linkpar.StreetWidth(1)); %same street width for all links

% WIM parameters, common to all links
SampleDensity=wimpar.SampleDensity;
NumTimeSamples=wimpar.NumTimeSamples;
%N=wimpar.NumPaths;
M=wimpar.NumSubPathsPerPath;
CenterFrequency=wimpar.CenterFrequency;
DelaySamplingInterval=wimpar.DelaySamplingInterval;
PathLossModel=wimpar.PathLossModel;
RandomSeed=wimpar.RandomSeed;
UniformTimeSampling=wimpar.UniformTimeSampling;
PathLossModelUsed=wimpar.PathLossModelUsed;
ShadowingModelUsed=wimpar.ShadowingModelUsed;
AnsiC_core=wimpar.AnsiC_core;
LookUpTable=wimpar.LookUpTable;
FixedPdpUsed = wimpar.FixedPdpUsed;
FixedAnglesUsed = wimpar.FixedAnglesUsed;
PolarisedArrays = wimpar.PolarisedArrays;

% antenna parameters
BsGainPattern=antpar.BsGainPattern;
BsGainAnglesAz=antpar.BsGainAnglesAz;
BsElementPosition=antpar.BsElementPosition;
MsGainPattern=antpar.MsGainPattern;
MsGainAnglesAz=antpar.MsGainAnglesAz;
MsElementPosition=antpar.MsElementPosition;
InterpFunction=antpar.InterpFunction;
InterpMethod=antpar.InterpMethod;

% link parameters
ScenarioVector = linkpar.ScenarioVector;
PropagConditionVector = linkpar.PropagConditionVector;
MsBsDistance=linkpar.MsBsDistance;
ThetaBs=linkpar.ThetaBs;
ThetaMs=linkpar.ThetaMs;
MsVelocity=linkpar.MsVelocity;
MsDirection=linkpar.MsDirection;
StreetWidth=linkpar.StreetWidth;
NumFloors = linkpar.NumFloors;
Dist1=linkpar.Dist1;


% extract the number of links
NumLinks=length(MsBsDistance);

% Check that the struct linkpar has the same number of parameters in
% each of its fields. This is also the number of links/users.
if (    NumLinks ~= length(ThetaBs)     ||...
        NumLinks ~= length(ThetaMs)     ||...
        NumLinks ~= length(MsVelocity)  ||...
        NumLinks ~= length(MsDirection)  ||...
        NumLinks ~= length(StreetWidth)  ||...
        NumLinks ~= length(NumFloors)  ||...
        NumLinks ~= length(Dist1))
    
    error('All fields in input struct LINKPAR must be of same size!')
end

% If layout parameters are defined, check for consistency
if (    isfield(linkpar,'BsXY')    ||...
        isfield(linkpar,'NofSect') ||...
        isfield(linkpar,'BsOmega') ||...
        isfield(linkpar,'MsXY')    ||...
        isfield(linkpar,'MsOmega') ||...
        isfield(linkpar,'Pairing')   )
    if (    size(linkpar.BsXY,2)~=size(linkpar.NofSect,2)    ||...
            size(linkpar.Pairing,1)~=sum(linkpar.NofSect)    ||...
            size(linkpar.Pairing,2)~=size(linkpar.MsXY,2)    ||...
            size(linkpar.BsOmega,1)~=max(linkpar.NofSect)    ||...
            size(linkpar.MsOmega,2)~=size(linkpar.MsXY,2)    ||...
            sum(linkpar.Pairing(:))~=NumLinks     )
        error('Layout parameters are inconsistent! See help LAYOUTPARSET.')
    end
end

% Set random seeds if given
if (isempty(RandomSeed)==0)
    rand('state',RandomSeed);
    randn('state',RandomSeed);
else
    rand('state',sum(100*clock));
    randn('state',sum(101*clock));
end



% determine the size of the MIMO system
% S - number of BS array antenna elements
if (numel(BsGainPattern)==1)
    S=wimpar.NumBsElements;
else
    S=size(BsGainPattern,1);
end

% U - number of MS array antenna elements
if (numel(MsGainPattern)==1)
    U=wimpar.NumMsElements;
else
    U=size(MsGainPattern,1);
end

% check that element displacement vector is of right size
if (length(BsElementPosition)~=S && length(BsElementPosition)~=1)
    error('antpar.BsElementPosition has wrong size!')
end

if (length(MsElementPosition)~=U && length(MsElementPosition)~=1)
    error('antpar.MsElementPosition has wrong size!')
end


% check that LUT size is a power-of-two
% this check is now also in wim_mex_core.c
if (strcmpi(AnsiC_core,'yes')==1)
    if (LookUpTable>0)
        if (2^nextpow2(LookUpTable)-LookUpTable~=0)
            wimpar.LookUpTable=2^nextpow2(LookUpTable);
            warning('MATLAB:LUTSizeChanged',['wimpar.LookUpTable is not a power-of-2: size changed to ' num2str(wimpar.LookUpTable) '.'])
        end
    end
end


if (strcmpi(wimpar.IntraClusterDsUsed,'yes')==1) & (strcmpi(AnsiC_core,'yes')==1)
    warning('AnsiC_core does not support the IntraClusterDs yet! wimpar.AnsiC_core set to "no"')
    wimpar.AnsiC_core = 'no';
end


% GENERATION OF RANDOM "BULK" PARAMETERS FOR ALL LINKS
switch (ni)

    case (3)    % do the basic thing

        % generate bulk parameters for all links
        %bulkpar=generate_bulk_par_polarised(wimpar,linkpar,antpar,fixpar);
        bulkpar=generate_bulk_par(wimpar,linkpar,antpar,fixpar);

        % get number of clusters from bulk parameters (located here because
        %  for the case of FixedPdpUsed
        N = size(bulkpar.delays,2);

        % for interpolation
        aods=bulkpar.aods;
        aoas=bulkpar.aoas;


    case (4)    % do not generate random link parameters, use initial values

        % take bulk parameters from input struct
        bulkpar=initvalues;
        
        % This IF is added to remove intra cluster delay spred effects from
        % initial values (spread takes effect in wim_core.m)
        if strcmp(wimpar.IntraClusterDsUsed,'yes')
            for k=1:NumLinks
                % Remove intra cluster delay values from initial values
                tmp = initvalues.delays(k,:);
                tmp([initvalues.IndexOfDividedClust(k,1)+[1:2],initvalues.IndexOfDividedClust(k,2)+2+[1:2]])= [];
                tmpDelay(k,:) = tmp;
                % Remove intra cluster power values from initial values
                tmp = initvalues.path_powers(k,:);
                tmp([initvalues.IndexOfDividedClust(k,1)+[1:2],initvalues.IndexOfDividedClust(k,2)+2+[1:2]])= [];
                tmpPower(k,:) = tmp;
            end
            bulkpar.delays = tmpDelay;
            bulkpar.path_powers = tmpPower;
        end

        % get number of clusters from bulk parameters (located here because
        %  for the case of FixedPdpUsed
        N = size(bulkpar.delays,2);

        % for interpolation
        aods=bulkpar.aods;
        aoas=bulkpar.aoas;


end



% ANTENNA FIELD PATTERN INTERPOLATION
% Interpolation is computationally intensive, so avoid it if possible.
% Since elevation will not be supported, dismiss the elevation dimension (for now)
% NOTE: aods/aoas should be given in degrees.
BsGainIsScalar=0;
MsGainIsScalar=0;
if numel(BsGainPattern)>1
    if strcmpi(PolarisedArrays,'yes')
        BsGainPatternInterpolated = zeros([2 S size(aods)]); % [polarizations(2) elements links N(6) M(20)]
        BsGainPatternInterpolated(1,:,:,:,:)=feval(InterpFunction,squeeze(BsGainPattern(:,1,1,:)),BsGainAnglesAz,aods, InterpMethod); % V
        BsGainPatternInterpolated(2,:,:,:,:)=feval(InterpFunction,squeeze(BsGainPattern(:,2,1,:)),BsGainAnglesAz,aods, InterpMethod); % H
        BsGainPatternInterpolated=permute(BsGainPatternInterpolated,[3 2 1 4 5]); % [link rx_element polarization path subpath]
    else
        BsGainPatternInterpolated=feval(InterpFunction,squeeze(BsGainPattern(:,1,1,:)),BsGainAnglesAz,aods, InterpMethod); % V only
        BsGainPatternInterpolated=permute(BsGainPatternInterpolated,[2 1 3 4]);
    end
else    % if BsGainPattern is scalar
    if strcmpi(PolarisedArrays,'yes')
        BsGainPatternInterpolated=repmat(BsGainPattern, [NumLinks S 2 N M]);    % [link rx_element polarization path subpath]
        BsGainIsScalar=1;
    else
        BsGainPatternInterpolated=repmat(BsGainPattern, [NumLinks S N M]);
        BsGainIsScalar=1;
    end
end

if numel(MsGainPattern)>1
    if strcmpi(PolarisedArrays,'yes')
        MsGainPatternInterpolated=zeros([2 U size(aoas)]);% [polarizations(2) elements links N(6) M(20)]
        MsGainPatternInterpolated(1,:,:,:,:)=feval(InterpFunction,squeeze(MsGainPattern(:,1,1,:)),MsGainAnglesAz,aoas, InterpMethod); % V
        MsGainPatternInterpolated(2,:,:,:,:)=feval(InterpFunction,squeeze(MsGainPattern(:,2,1,:)),MsGainAnglesAz,aoas, InterpMethod); % H
        MsGainPatternInterpolated=permute(MsGainPatternInterpolated,[3 2 1 4 5]); % [link Ms_element polarization path subpath]
    else
        MsGainPatternInterpolated=feval(InterpFunction,squeeze(MsGainPattern(:,1,1,:)),MsGainAnglesAz,aoas, InterpMethod); % V only
        MsGainPatternInterpolated=permute(MsGainPatternInterpolated,[2 1 3 4]);
    end
else    % if MsGainPattern is scalar
    if strcmpi(PolarisedArrays,'yes')
        MsGainPatternInterpolated=repmat(MsGainPattern, [NumLinks U 2 N M]);    % [link rx_element polarization path subpath]
        MsGainIsScalar=1;
    else
        MsGainPatternInterpolated=repmat(MsGainPattern, [NumLinks U N M]);
        MsGainIsScalar=1;
    end
end

% Note: The gain patterns at this point have size(MsGainPatternInterpolated) = [link rx_element path subpath]
%  OR size(MsGainPatternInterpolated) = [link rx_element polarization path subpath] (the same for BsGainPatternInterpolated)


%% Do antenna field pattern interpolation for the LOS path
%  Note! this is done even for NLOS (but result is not used)
% Polarised arrays case added 19.12.2005, PekKy

% BS antenna
if numel(BsGainPattern)>1
    BsGain_Theta_BS= feval(InterpFunction,squeeze(BsGainPattern(:,1,1,:)),BsGainAnglesAz,ThetaBs(:), InterpMethod); % V only
    BsGain_Theta_BS= BsGain_Theta_BS.'; % size()= [NumLinks S]
    %         if strcmpi(PolarisedArrays,'yes')
    %             tmp = feval(InterpFunction,squeeze(BsGainPattern(:,2,1,:)),BsGainAnglesAz,ThetaBs(:), InterpMethod); % H pol
    %             BsGain_Theta_BS(:,:,2) = tmp.'; % size()= [NumLinks S 2]
    %         end
else
    BsGain_Theta_BS=repmat(BsGainPattern,[NumLinks S]);
    %         if strcmpi(PolarisedArrays,'yes')
    %             BsGain_Theta_BS(:,:,2) = BsGain_Theta_BS;   % H pol
    %         end
end

% MS antenna
if numel(MsGainPattern)>1
    MsGain_Theta_MS= feval(InterpFunction,squeeze(MsGainPattern(:,1,1,:)),MsGainAnglesAz,ThetaMs(:), InterpMethod); % V only
    MsGain_Theta_MS= MsGain_Theta_MS.'; % size()= [NumLinks U]
    %         if strcmpi(PolarisedArrays,'yes')
    %             tmp = feval(InterpFunction,squeeze(MsGainPattern(:,2,1,:)),MsGainAnglesAz,ThetaMs(:), InterpMethod); % H pol
    %             MsGain_Theta_MS(:,:,2) = tmp.'; % size()= [NumLinks U 2]
    %         end
else
    MsGain_Theta_MS= repmat(MsGainPattern,[NumLinks U]);
    %         if strcmpi(PolarisedArrays,'yes')
    %             MsGain_Theta_MS(:,:,2) = MsGain_Theta_MS;   % H pol
    %         end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Channel Matrix Generation    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Separate processing for LOS and NLOS scenarios %%%%
%
if sum(bulkpar.propag_condition)>0      % LOS links

    PCind = find(bulkpar.propag_condition);   % Propagation condition index, LOS

    % CHANNEL MATRIX GENERATION
    [Htmp delta_t FinalPhases FinalPhases_LOS] = wim_core( wimpar,...
        linkpar,...
        antpar,...
        bulkpar,...
        BsGainPatternInterpolated,...
        BsGain_Theta_BS,...             % gain of LOS path
        MsGainPatternInterpolated,...
        MsGain_Theta_MS,...             % gain of LOS path
        0,...                           % offset time (not used typically)
        BsGainIsScalar,...
        MsGainIsScalar,...
        PCind);

    H(:,:,:,:,PCind) = Htmp;

    % final phases
    bulkpar.subpath_phases(PCind,:,:,:)=FinalPhases;
    % time sampling grid
    bulkpar.delta_t(PCind)=delta_t;
end

if (length(bulkpar.propag_condition)-sum(bulkpar.propag_condition))>0   % NLOS links

    PCind = find(bulkpar.propag_condition==0);   % Propagation condition index, NLOS

    % CHANNEL MATRIX GENERATION
    [Htmp delta_t FinalPhases] = wim_core( wimpar,...
        linkpar,...
        antpar,...
        bulkpar,...
        BsGainPatternInterpolated,...
        BsGain_Theta_BS,...             % gain of LOS path
        MsGainPatternInterpolated,...
        MsGain_Theta_MS,...             % gain of LOS path
        0,...                           % offset time (not used typically)
        BsGainIsScalar,...
        MsGainIsScalar,...
        PCind);

    H(:,:,:,:,PCind) = Htmp;

    % final phases
    bulkpar.subpath_phases(PCind,:,:,:)=FinalPhases;
    % time sampling grid
    bulkpar.delta_t(PCind)=delta_t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%


% If path loss and shadowing are to be multiplied into the output
if ( strcmpi(PathLossModelUsed,'yes') || strcmp(ShadowingModelUsed,'yes') )

    O2VLinks = [or(~isnan(bulkpar.LoS02VLinks), ~isnan(bulkpar.NLoS02VLinks))]';
    O2Vshadowing = bulkpar.o2v_shadow_fading.*O2VLinks;

    if (size(H,5)==1) % only one link
        if strcmpi(PathLossModelUsed,'yes')
            H=sqrt(bulkpar.path_losses).*H;   % path loss in linear scale
        end

        if strcmpi(ShadowingModelUsed,'yes')
            H=H*sqrt(bulkpar.shadow_fading)*sqrt(O2Vshadowing);   % shadow fading in linear scale
        end
    else    % if more than one link

        siz_H=size(H);
        Hmat=reshape(H,prod(siz_H(1:end-1)),siz_H(end));  % a matrix with NumLinks cols
        if strcmpi(PathLossModelUsed,'yes')
            pl_mat=diag(sparse(sqrt(bulkpar.path_losses)));
            Hmat=Hmat*pl_mat;           % multiply path loss into each link
        end


        if strcmp(ShadowingModelUsed,'yes')
            sf_mat = diag(sparse(sqrt(bulkpar.shadow_fading)));    % shadow fading is in linear scale
            sf_mat2 = diag(sparse(sqrt(O2Vshadowing)));    % shadow fading is in linear scale
            Hmat=Hmat*sf_mat*sf_mat2;         % multiply shadow fading into each link
        end

        H=reshape(Hmat,siz_H);      % put back to original size

    end
end


% GENERATE OUTPUT
no=nargout;

if strcmpi(wimpar.IntraClusterDsUsed,'yes')
    bulks = bulkpar;
    bulkpar.delays = repmat(NaN,size(bulks.delays,1),size(bulks.delays,2)+4);
    bulkpar.path_powers = repmat(NaN,size(bulks.delays,1),size(bulks.delays,2)+4);
    for link = 1:NumLinks
        B5ind = find((linkpar.ScenarioVector>=7 & linkpar.ScenarioVector<=9));
        if link==B5ind
            bulkpar.delays(link,1:length(bulks.delays(link,:))) = bulks.delays(link,:);
            bulkpar.path_powers(link,1:length(bulks.delays(link,:))) = bulks.path_powers(link,:);
        else
            P = bulks.path_powers(link,:); P(isnan(P)) = -Inf;
            SortedPower = fliplr(sort(P,2)); P(isinf(P)) = NaN;
            %SubClustInd = P > SortedPower(3); % Index of the cluster to be divided
            % Find index to the clusters to be divided
            [tmp tmpind] = sort(P); SubClustInd = tmpind(end-1:end);
            SubClustInd = zeros(1,size(P,2)); SubClustInd(tmpind(end-1:end)) = 1;
            SubClustDelays = [0  5e-009  10e-009]';
            taus = repmat(bulks.delays(link,:),length(SubClustDelays),1) + repmat(SubClustDelays,1,length(P));
            taus(2:3,~SubClustInd) = NaN;
            taus = reshape(taus,1,length(SubClustDelays)*size(P,2));
            taus(isnan(taus)) = [];
            powers = repmat(bulks.path_powers(link,:),length(SubClustDelays),1);
            SubClustP = [10/20 6/20 4/20]';
            powers(:,find(SubClustInd==1)) = powers(:,find(SubClustInd==1)).*repmat(SubClustP,1,2);
            powers(2:3,~SubClustInd) = NaN;
            powers = reshape(powers,1,length(SubClustP)*size(P,2));
            powers(isnan(powers)) = [];
            IndexOfDividedClust(link,:) = find(SubClustInd==1);


            if (no>1)
                delays(link,1:length(taus)) = taus;
                if (no>2)
                    bulkpar.delays(link,1:length(taus)) = taus;
                    bulkpar.path_powers(link,1:length(taus)) = powers;
                    bulkpar.IndexOfDividedClust = IndexOfDividedClust;
                    bulkpar.aods = aods;
                    bulkpar.aoas = aoas;
                    full_output = bulkpar;
                end
            end
            clear taus powers
        end
    end

else
    if (no>1)
        delays = bulkpar.delays;
        if (no>2)
            if (sum(bulkpar.propag_condition)>0)     % At least one LOS links included
                bulkpar.Phi_LOS=FinalPhases_LOS;
                full_output=bulkpar;
            else    % Only NLOS links included
                full_output=bulkpar;
            end
        end
    end
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);