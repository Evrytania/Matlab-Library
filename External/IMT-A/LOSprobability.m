function indLOS=LOSprobability(wimpar,linkpar,fixpar,iterpar)
%LOSPROBABILITY Random LOS/NLOS condition generation for WIM
%   INDLOS=LOSPROBABILITY(WIMPAR,LINKPAR) is a vector defining links
%   propagation condition. Vector elements have values '0' or '1' and length is
%   number of links. '0' stands for NLOS and '1' stands for LOS link,
%   vector length is number of links. LOS/NLOS condition is drawn randomly
%   for each link according to LOS probabilities defined in [1, Table 4-7].
%
%   Ref. [1]: D1.1.1 V1.0, "WINNER II interim channel models"
%        [2]: IMT.EVAL
%
%   See also WIM, WIMPARSET, LINKPARSET, ANTPARSET.

%   Authors: Pekka Ky�sti (EBIT), Mikko Alatossava (CWC/UOULU)

%   Modifications: Parameters from [1, Table 4-7] are used:         15.1.2007 Marko
%                  Divider for C1 changed to 200, used to be 500    12.2.2007 MikkoA
%                  Updated based on [2]                             22.9.2007 PekKy
%                  Update InH to match with Seoul contribution      25.9.2008 Zhiwen Wu  

NumLinks = length(iterpar.UserIndeces);
Scenario = iterpar.Scenario;
MsBsDistance = linkpar.MsBsDistance(iterpar.UserIndeces);

% Probability of LOS
switch upper(Scenario)      % See Table A1-3 in [2]
    
    case {'A2'}     % InH
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 3*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        % Distance between MS and BS locations (on ground level)
        r = sqrt(MsBsDistance.^2-(BsHeight-MsHeight).^2);
        % Calculate LOS probability
        pLOS = exp(-(r-18)/27);
        pLOS(r<=18) = 1;         % If r<=18,  pLOS = 1
        pLOS(r>=37) = 0.5;        % If r>=37, pLOS = 0.5
                
    case {'B1'}     % UMi
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 10*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        % Distance between MS and BS locations (on ground level)
        r = sqrt(MsBsDistance.^2-(BsHeight-MsHeight).^2);
        % Calculate LOS probability
        pLOS = min(18./r,ones(size(r))).*(1-exp(-r/36)) + exp(-r/36);

    case {'C1'}     % SMa
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 35*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        % Distance between MS and BS locations (on ground level)
        r = sqrt(MsBsDistance.^2-(BsHeight-MsHeight).^2);
        pLOS = exp(-(r-10)/200);
        pLOS(r<10) = 1;     % If r<10, pLOS = 1
        
    case {'C2'}     % UMa
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 25*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        % Distance between MS and BS locations (on ground level)
        r = sqrt(MsBsDistance.^2-(BsHeight-MsHeight).^2);
        pLOS = min(18./r,ones(size(r))).*(1-exp(-r/63)) + exp(-r/63);
        
    case {'D1'}     % RMa
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 35*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        % Distance between MS and BS locations (on ground level)
        r = sqrt(MsBsDistance.^2-(BsHeight-MsHeight).^2);
        pLOS = exp(-(r-10)/1000);
        pLOS(r<10) = 1;     % If r<10, pLOS = 1
        
        
end

% output, 0 for NLOS and 1 for LOS links
indLOS = rand(1,NumLinks)<pLOS;
