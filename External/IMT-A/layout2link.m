function linkpar=layout2link(layoutpar)
%LAYOUT2LINK Layout to link parameter conversion for WIM.
%   LINKPAR=LAYOUT2LINK(LAYOUTPAR) returns extended set of link parameters
%   in the case of layout parameters are defined. It converts layout
%   parameters to Ms/Bs distances, LOS directions etc. LAYOUT2LINK is used with 
%   WIM the following way: [..] = wim(wimparset,layout2link(layoutpar),antpar).
%   The function pair LAYOUT2LINK(LAYOUTPAR(..)) is optional to function
%   LINKPARSET and will be used if MS/BS locations are given in (x,y)
%   co-ordinates.
%
%   The output parameters are:
%
%   BsXY         - directly from LAYOUTPAR, NOTE! rounding to 1m resolution 
%   NofSect      - directly from LAYOUTPAR 
%   BsOmega      - directly from LAYOUTPAR
%   MsXY         - directly from LAYOUTPAR, NOTE! rounding to 1m resolution 
%   MsOmega      - directly from LAYOUTPAR
%   Pairing      - directly from LAYOUTPAR
%   MsBsDistance - calculated from LAYOUTPAR
%   ThetaBs      - calculated from LAYOUTPAR
%   ThetaMs      - calculated from LAYOUTPAR
%   MsVelocity   - defined in LAYOUTPAR 
%   MsDirection  - defined in LAYOUTPAR
%   StreetWidth  - defined in LAYOUTPAR
%   Dist2        - defined in LAYOUTPAR
%
%   See [1, Fig. 6.1 and 6.2].
%
%   Ref. [1]: D1.1.1 V1.0, "WINNER II interim channel models"
%
%   See also WIM, LAYOUTPARSET, WIMPARSET, LINKPARSET, ANTPARSET.

%   Authors: Pekka Kyösti (EBIT)
%


BsXY = round(layoutpar.BsXY);
NofSect = layoutpar.NofSect;
BsOmega = layoutpar.BsOmega;
MsXY = round(layoutpar.MsXY);
MsOmega = layoutpar.MsOmega;
Pairing = layoutpar.Pairing;
 
% linkpar struct with layout parameters included
linkpar=struct( 'BsXY',    [],...
                'NofSect', [],...
                'BsOmega', [],...
                'MsXY',    [],...
                'MsOmega', [],...   
                'Pairing', [],...
                'ScenarioVector',[],...
                'PropagConditionVector',[],...
                'MsBsDistance', [],...
                'ThetaBs', [],...
                'ThetaMs', [],...
                'MsHeight',[],... 
                'BsHeight',[],...
                'MsVelocity', [],...
                'MsDirection', [],...
                'StreetWidth', [],...
                'NumFloors', [],...
                'NumPenetratedFloors', [],...
                'Dist1', []);


% get number of links and Ms & Sector indices            
K = sum(Pairing(:));    % number of links
[r c] = find(Pairing);
indMs = c';      % links coupling to MSs
indSect = r';    % links coupling to Sectors
% get links coupling to BSs
tmpsum = cumsum(NofSect);
for k=1:length(r)
    tmp = find(diff(r(k)<=tmpsum));
    if isempty(tmp)
        indBs(k)=1;
    else
        indBs(k)=tmp+1;
    end
end


% MS-Sector distance
linkpar.MsBsDistance = sqrt(sum((BsXY(:,indBs)-MsXY(:,indMs)).^2));  % [1, eq. 6.2]

% LOS direction from BS Sector array to MS wrt array broad side
SectOmega = BsOmega(indSect);
if size(SectOmega,1)>1 SectOmega = SectOmega'; end
%tmp = -atan((MsXY(2,indMs)-BsXY(2,indBs))./(MsXY(1,indMs)-BsXY(1,indBs)))-SectOmega;
tmp = -atan((MsXY(2,indMs)-BsXY(2,indBs))./(MsXY(1,indMs)-BsXY(1,indBs)))/pi*180-SectOmega; 
ThetaBs = tmp+((MsXY(1,indMs)>=BsXY(1,indBs))*2-1)*90;   % [1, eq. 6.2]
linkpar.ThetaBs = prin_value(ThetaBs);

% LOS direction from MS array to BS wrt array broad side
%tmp = -atan((BsXY(2,indBs)-MsXY(2,indMs))./(BsXY(1,indBs)-MsXY(1,indMs)))-MsOmega(indMs);
tmp = -atan((BsXY(2,indBs)-MsXY(2,indMs))./(BsXY(1,indBs)-MsXY(1,indMs)))/pi*180-MsOmega(indMs);
ThetaMs = tmp+((BsXY(1,indBs)>=MsXY(1,indMs))*2-1)*90;   % similar to [1, eq. 6.2]
linkpar.ThetaMs = prin_value(ThetaMs);

% the rest of the link parameters
linkpar.ScenarioVector        = layoutpar.ScenarioVector;
linkpar.PropagConditionVector = layoutpar.PropagConditionVector;
linkpar.MsHeight              = layoutpar.MsHeight(indMs);
linkpar.BsHeight              = layoutpar.BsHeight(indBs);
linkpar.MsVelocity            = layoutpar.MsVelocity(indMs);
linkpar.MsDirection           = layoutpar.MsDirection(indMs);
linkpar.StreetWidth           = layoutpar.StreetWidth;
linkpar.NumFloors             = layoutpar.NumFloors;
linkpar.NumPenetratedFloors   = layoutpar.NumPenetratedFloors;
linkpar.Dist1                 = layoutpar.Dist1;

% convert actual layoutpar to linkpar
linkpar.BsXY    = BsXY;
linkpar.NofSect = NofSect;
linkpar.BsOmega = BsOmega;
linkpar.MsXY    = MsXY;
linkpar.MsOmega = MsOmega;
linkpar.Pairing = Pairing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);
