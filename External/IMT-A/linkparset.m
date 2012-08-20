function linkpar=linkparset(varargin)
%LINKPARSET Link parameter configuration for WIM
%   LINKPAR=LINKPARSET(K) is a struct consisting of randomly generated link
%   parameters for K links. LINKPAR=LINKPARSET(K,RMAX) uses cell radius
%   RMAX for generation of MS-BS distances (default: 100 meters).
%
%   LINKPAR=LINKPARSET(...,SEED) sets the random seed used in link
%   parameter generation. 
%
%   The parameters and their defaults are:
%
%   ScenarioVector        - maps scenario names to links (see below)
%   PropagConditionVector - maps propagation condition (NLOS=0/LOS=1) to links
%   MsBsDistance          - see below
%   BsHeight              - NaN default (see below)
%   MsHeight              - NaN default (see below)
%   ThetaBs               - U(-180,180) degrees, U denotes uniform pdf
%   ThetaMs               - U(-180,180) degrees
%   MsVelocity            - 10 meters per second 
%   MsDirection           - U(-180,180) degrees with respect to broadside
%   StreetWidth           - Average width of the streets, same for all users.
%   LayoutType            - Layout type for UMi (B1/B4) path loss, {0}=hexagonal, 1=Manhattan
%   NumFloors             - For scenarios A2/B4 this determines the floor number where BS/MS is. 
%   OtoI_OutdoorPL        - Outdoor-to-Indoor propagation condition for UMi (see. note 3 in A1-2, [1])
%   Dist1                 - Distance from BS to "the street crossing" (last LOS point), default NaN -> will be drawn randomly
%   BuildingHeight        - Average building height, same for all users

%  
%   The pdf of the random variable (RV) MsBsDistance is R+RMIN, where R is
%   an RV with pdf p(r)=2*r/r0^2, where r0 defaults to (RMAX-RMIN) meters.
%   Hence, MsBsDistance is an RV such that users are approximately
%   uniformly distributed in a circular disk over [RMIN,RMAX] meters. Default
%   values for RMIN and RMAX are 5 and 100 m respectively. For usability,
%   the same default is used for all scenarios. 
%  NOTE! IF random MsBsDistance is used, SET FEASIBLE values for
%   RMIN and RMAX. This is not done automatically for different scenarios.
%   
%   Some notes about the parameters:
%
%   - Bs/MsHeight, if(isnan) -> default heights are used from [3,table 4-4].  
%   - StreetWidth, this is utilized only with path loss model in [1, sec 5.4.1.2]
%   - Dist1 is defined in [1, Figure 4-3] and generated randomly if empty
%     
%   Ref. [1]: D1.1.2 V1.0, "WINNER II channel models"
%        [2]: 3GPP TR 25.996 v6.1.0 (2003-09)
%        [3]: D1.1.1 V1.2, "WINNER II interim channel models"
%        [4]: IMT-EVAL.
%
%   See also WIM, WIMPARSET, LAYOUTPARSET, ANTPARSET.

%   Authors: Jari Salo (HUT), Pekka Kyösti (EBIT), Daniela Laselva (EBIT), 
%   Giovanni Del Galdo (TUI), Marko Milojevic (TUI), Christian Schneider (TUI)
%   Lassi Hentilä (EBIT), Mikko Alatossava (CWC/UOULU)


% defaults
num=1;      % number of links
rmax=100;   % cell radius [m]
rmin=50;     % minimum distance [m]

ni=length(varargin);
if ni>0, if (~isempty(varargin{1})), num=varargin{1}; end, end
if ni>1, if (~isempty(varargin{2})), rmax=varargin{2}; end, end
if ni>2, if (~isempty(varargin{3})), seed=varargin{3}; rand('state',floor(seed)); end, end
if ni>3, error('Too many input arguments!'), end

%Scenario need to be given for each user in (1 x K) size vector in numerical format.
%From this vector the scenario for each user is mapped as 
% 1=InH (A2), 2=UMi (B1), 34=SMa (C1), 4=UMa (C2), 5=UMi O-2-I (B4), 6=RMa (D1)

            
linkpar=struct( 'ScenarioVector',1*ones(1,num),... % A2, B1, B4, C1, C2, D1
                'PropagConditionVector',zeros(1,num),...
                'MsBsDistance',distrnd(num,rmax-rmin)+rmin,...
                'BsHeight',repmat(NaN,1,num),... % NaN -> default heights from D1.1.2
                'MsHeight',repmat(NaN,1,num),... % NaN -> default heights from D1.1.2
                'ThetaBs',360*(rand(1,num)-0.5),...
                'ThetaMs',360*(rand(1,num)-0.5),...
                'MsVelocity',repmat(10,1,num),...
                'MsDirection',360*(rand(1,num)-0.5),...
                'StreetWidth',20*ones(1,num),...
                'LayoutType',0*ones(1,num),...  % Layout type for UMi (B1/C4) path loss, 0=hexagonal, 1=Manhattan
                'NumFloors',1*ones(1,num),... The ground floor is number 1
                'OtoI_OutdoorPL',1*ones(1,num),... Outdoor-to-Indoor propagation condition for UMi, 1 is LOS, 0 is NLOS (see. note 3 in A1-2, [4])
                'Dist1',repmat(NaN,1,num),...
                'BuildingHeight', repmat(NaN,1,num),... % NaN -> default heights from ScenParTables
                'LoS02ILinks', repmat(NaN,1,num),... % NaN -> will be drawn randomly
                'LoS02VLinks', repmat(NaN,1,num),... % NaN -> will be drawn randomly
                'NLoS02ILinks', repmat(NaN,1,num),... % NaN -> will be drawn randomly
                'NLoS02VLinks', repmat(NaN,1,num)); % NaN -> will be drawn randomly
            
            
function d=distrnd(num,rmax)
% DISTRND Distance from BS in a circular cell
%   D=DISTRND(K,RMAX) generates K random variables from the pdf
%   p(r)=2*r/RMAX^2. This is the pdf of distance from base station when
%   users are uniformly (in area) distributed in a cell with radius RMAX.

%   Authors: Jari Salo (HUT), Marko Milojevic (TUI) 
%   $Revision: 0.1$  $Date: Sep 30, 2004$

% create random variables from triangular pdf whose width is 2*rmax
a=sum(rmax*rand(2,num));

% fold the random variables about the rmax
inds=find(a>rmax);
a(inds)=-a(inds)+2*rmax;

d=a(:).';