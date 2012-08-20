function layoutpar=layoutparset(varargin)
%LAYOUTPARSET Link parameter configuration for WIM
%   LAYOUTPAR=LAYOUTPARSET(NofMs,NofBs,SectPerBs,K) is a struct consisting 
%   of randomly generated network layout parameters. BS and MS positions 
%   are set and a pairing matrix with K links between is generated.
%   LAYOUTPAR=LAYOUTPARSET(NofMs,NofBs,K,RMAX) uses layout range
%   RMAX for generation of MS and BS positions on cartesian co-ordinate
%   system (default: 100 meters).
%
%   LAYOUTPAR=LAYOUTPARSET(...,SEED) sets the random seed used in layout
%   parameter generation. 
%
%   The parameters and their defaults are:
%
%   BsXY        - matrix of BS (x,y) co-ordinates with dimensions 2xNofBs [m]
%   NofSect     - vector of number of sectors in each of the BSs, default=ones(1,NofBs)
%   BsOmega     - matrix of BS array broad side orientations in a sector [deg]
%   MsXY        - matrix of MS (x,y) co-ordinates with dimensions 2xNofMs  [m]
%   MsOmega     - vector of MS array broad side orientations [deg] [1, Fig 5-2]
%   Pairing     - matrix defining which links are modelled, NofSect x NofMs
%
%   ScenarioVector        - maps scenario names to links (see ScenarioMapping.m)
%   PropagConditionVector - maps propagation condition (NLOS=0/LOS=1) to links
%   MsVelocity  - 10 meters per second 
%   MsDirection - U(-180,180) degrees with respect to broadside
%   StreetWidth - 25 meters
%   NumFloors   - For scenarios A2/B4 this determines the floor number of BS/MS 
%   NumPenetratedFloors   - Number of floor between BS/MS for the A1 path loss (default is zero)
%   Dist1       - Distance from BS to "the street crossing" (last LOS point), default NaN -> will be drawn randomly

%
%   See [1, Fig. 5-1 and 5-2].
%   
%   Some notes about the parameters:
%
%   - Co-ordinates of Bs and Ms should be given in meters with resolution
%     of 1 meter. One meter resolution is assumed in auto-correlation
%     generation of LScorrelation.m function.
%   - BsOmega is a matrix with dimensions max(NofSect)xNofBs. Each column
%     of the matrix contains orientations of sectorised arrays with respect
%     to some fixed North direction. If some BS have less sectors than
%     others, the non-existing sector orientations are set to zero. E.g.
%     setup with 2 BS, one with 1 sector and other with 3 sectors. In this
%     case orientation matrix could be e.g. BsOmega = [11 22; 0 33; 0 44].
%   - Pairing is a matrix with dimensions NofSect x NofMs, i.e. one entry for
%     each BS sector/MS pair. Value '1' stands for "link will be modelled"
%     and value '0' stands for "link will not be modelled". E.g. with all ones
%     matrix, all the MS are connected to all sectors. With e.g. first rows
%     ones and others zeros means, that all MS are connected to only 1st sector
%     of 1st BS.
%   - StreetWidth, this is utilized only with path loss model in [1, Table 4-4]
%   - Dist1 is defined in [1, Figure 4-3] and generated randomly if empty
%
%   Ref. [1]: D1.1.2 V1.0, "WINNER II channel models"
%
%   See also WIM, LAYOUT2LINK, WIMPARSET, LINKPARSET, ANTPARSET.

%   Authors: Pekka Kyösti (EBIT)
%
%   Updates to Phase II model:  Added new (linkpar) parameters ScenarioVector,
%              PropagConditionVector, NumFloors             (29.5.07  PekKy)


% defaults
NofMs=1;       % number of BMs
NofBs=1;       % number of BSs
K=1;           % number of links
rmax=500;      % layout range [m]
SectPerBs = 1; % default number of sectors in a BS
%BSrmin=10;     % minimum distance of BSs [m]

% inputs
ni=length(varargin);
if ni>0, if (~isempty(varargin{1})), NofMs=varargin{1}; end, end
if ni>1, if (~isempty(varargin{2})), NofBs=varargin{2}; end, end
if ni>2, if (~isempty(varargin{3})), SectPerBs=varargin{3}; end, end
if ni>3, if (~isempty(varargin{3})), K=varargin{4}; end, end
if ni>4, if (~isempty(varargin{4})), rmax=varargin{5}; end, end
if ni>5, if (~isempty(varargin{5})), seed=varargin{6}; rand('state',floor(seed)); end, end
if ni>6, error('Too many input arguments!'), end

% check input SectPerBs
if length(SectPerBs)==1 NofSect=repmat(SectPerBs,1,NofBs); % SectPerBs [scalar]: NofSect equal in any BS
elseif length(SectPerBs)~=NofBs                            % SectPerBs missing in some BSs
    SectPerBs = 1; NofSect=repmat(SectPerBs,1,NofBs);        % set default values
    warning('Number of sectors required for each Bs')
    disp(['SectPerBs set to default ' mat2str(SectPerBs)])
else                                                       % SectPerBs [1*NofBs], NofSect differ among BSs
    NofSect=SectPerBs;
end

% calculate output BsOmega 
BsOmega=[];
for i=1:NofBs 
    alpha = 360*rand;            % random broadside angle of the first Bs antenna
    OmegaTmp=(360/NofSect(i))*repmat([1:NofSect(i)]',1)+alpha; % BsOmega of a Bs
    if length(OmegaTmp)<max(NofSect) for tmp=length(OmegaTmp)+1:max(NofSect) OmegaTmp=[OmegaTmp; 0]; end,end
    BsOmega=[BsOmega OmegaTmp]; 
    BsOmega=prin_value(BsOmega); % BsOmega related to the North Reference
end 

% outputs
layoutpar=struct('BsXY',    round(rand(2,NofBs)*rmax),...
                 'NofSect', NofSect,...
                 'BsOmega', BsOmega,...
                 'MsXY',    round(rand(2,NofMs)*rmax),...
                 'MsOmega', 360*(rand(1,NofMs)-0.5),...   
                 'Pairing', fillpairing(NofMs,NofBs,SectPerBs,K),...
                 'ScenarioVector',1*ones(1,K),... % InH, UMi, SMa, SMa, UMiO2I, RMa
                 'PropagConditionVector',round(rand(1,K)),...
                 'MsHeight',1.5*ones(1,NofMs),... 
                 'BsHeight',32*ones(1,NofBs),... 
                 'MsVelocity',repmat(10,1,K),...        % linkparameters below this
                 'MsDirection',360*(rand(1,K)-0.5),...              
                 'StreetWidth',20*ones(1,K),...
                 'NumFloors', 1*ones(1,K),...   % The ground floor is number 1
                 'NumPenetratedFloors',0*ones(1,K),... Number of floor for the A1 path loss (default is zero)
                 'Dist1',repmat(NaN,1,K)); 
                
                

function A=fillpairing(NofMs,NofBs,NofSect,K)
% FILLPAIRING
%   A=FILLPAIRING(NOFMS,NOFBS,NOFSECT) generates pairing matrix A with
%   one link from each MS to a random sector. 

%   Authors: Pekka Kyösti (EBIT), Daniela Laselva (EBIT)

A=zeros(NofBs*max(NofSect),NofMs);
% generate '1' to K random entry of a matrix
if K>NofBs*max(NofSect)*NofMs
    K=NofBs*max(NofSect)*NofMs;
    warning('Number of modelled links limited by the layout')
    disp(['Number of links set to ' mat2str(K)])
end
tmp=randperm(NofBs*max(NofSect)*NofMs);
A(tmp(1:K))=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);
   