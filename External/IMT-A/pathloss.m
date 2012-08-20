function [loss, linkpar, fixpar, iterpar] = pathloss(wimpar,linkpar,fixpar,iterpar)
%PATHLOSS WIM pathloss models
%   PATH_LOSSES=PATHLOSS(WIMPAR,LINKPAR,FIXPAR,ITERPAR) returns path losses in dB scale
%   for all links defined in WIM input struct LINKPAR for the center
%   frequency and scenario given in WIMPAR. The output is a column vector
%   whose length is equal to the number of links defined in LINKPAR, e.g.
%   LENGTH(LINKPAR.MsBsDistance). The center frequencies and distances
%   must be specified in Herzes and meters, respectively.
%
%   Refs.   [1]: 3GPP TR 25.996 v6.1.0 (2003-09)
%           [2]: ITU-R IMT.EVAL Channel model, CF version May 2008
%
%   See also WIMPARSET, LINKPARSET and SCENPARTABLES.

%   Authors: Lassi Hentila (EBIT), Daniela Laselva (EBIT), Jari Salo (HUT),
%   Pekka Kyosti (EBIT), Marko Milojevic (TUI), Mikko Alatossava (CWC/UOULU)

%   Revision history after WIM release:
%   Frequency and antenna height dependence added.    4.5.2006 HentLas
%   B5 scenarios added.                               4.5.2006 HentLas
%   B1 LOS&NLOS updated.                              8.9.2006 PekKy
%   D1.1.1 parameters, new scenarios added            12.1.2007 Marko
%   LoS/NloS condition renewed, D1.1.1 parameters     8.2.2007 MikkoA
%   Updated according to the D1.1.2                   8.10.2007 HentLas
%   Parameters updated to comply with ITU-R IMT.EVAL  27.5.2008 PekKy
%   Updated according to the clarification doc.       14.4.2009 HentLas


% extract required parameters from the input structs
NumLinks = length(iterpar.UserIndeces);
MsBsDistance=linkpar.MsBsDistance(iterpar.UserIndeces);
Scenario=iterpar.Scenario;
NumFloors = linkpar.NumFloors(iterpar.UserIndeces);
CenterFrequency = wimpar.CenterFrequency*1e-9; % Center frequency -> GHz

PropagCondition = iterpar.PropagCondition;
LoSConnectionLinks = find(PropagCondition);
NumLoSConnectionLinks = length(LoSConnectionLinks);
NLoSConnectionLinks = find(PropagCondition==0);
NumNLoSConnectionLinks = length(NLoSConnectionLinks);
O2ILinks = round(rand(1,NumLinks)); % 50 percent of the users are indoors
O2VLinks = ~O2ILinks; % the rest 50 percent of the users are in a car
LoS02ILinks = find(and(O2ILinks,PropagCondition));
LoS02VLinks = find(and(O2VLinks,PropagCondition));
NLoSO2ILinks = find(and(O2ILinks,~PropagCondition));
NLoSO2VLinks = find(and(O2VLinks,~PropagCondition));

SF_sigma = [];
if ~isempty(LoSConnectionLinks)
    if MsBsDistance(LoSConnectionLinks) > iterpar.LoS.PL_range(2)
        error('MsBsDistance exceeds the maximum allowed cell radius, see IMT.EVAL Table 1')
    elseif MsBsDistance(LoSConnectionLinks) < iterpar.LoS.PL_range(1)
        error('MsBsDistance is below the minimum allowed cell radius, see IMT.EVAL Table 1')
    end
end

% if ~isempty(NLoSConnectionLinks)
%     if linkpar.LayoutType  % If 'Manhattan'
%         if MsBsDistance(NLoSConnectionLinks) > iterpar.NLoS.PL_range(2)
%             error('MsBsDistance exceeds the maximum allowed cell radius, see IMT.EVAL Table 1')
%         elseif MsBsDistance(NLoSConnectionLinks) < iterpar.NLoS.PL_range(1)
%             error('MsBsDistance is below the minimum allowed cell radius, see IMT.EVAL Table 1')
%         end
%     else % Hexagonal layout
%         if MsBsDistance(NLoSConnectionLinks) > iterpar.NLoS.PL_range_hex(2)
%             error('MsBsDistance exceeds the maximum allowed cell radius, see IMT.EVAL Table 1')
%         elseif MsBsDistance(NLoSConnectionLinks) < iterpar.NLoS.PL_range_hex(1)
%             error('MsBsDistance is below the minimum allowed cell radius, see IMT.EVAL Table 1')
%         end
%     end
% end

if ~isempty(NLoSConnectionLinks)
    if MsBsDistance(NLoSConnectionLinks) > iterpar.NLoS.PL_range(2)
        error('MsBsDistance exceeds the maximum allowed cell radius, see IMT.EVAL Table 1')
    elseif MsBsDistance(NLoSConnectionLinks) < iterpar.NLoS.PL_range(1)
        error('MsBsDistance is below the minimum allowed cell radius, see IMT.EVAL Table 1')
    end
    if linkpar.LayoutType & strcmpi(iterpar.Scenario,'B1')  % If 'Manhattan' and UMi (B1)
        if MsBsDistance(NLoSConnectionLinks) > iterpar.NLoS.PL_range_Manhattan
            error('MsBsDistance exceeds the maximum allowed cell radius, see IMT.EVAL Table 1')
        end
    end
end


if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
    warning('Bs/MsHeights not defined by the user --> defauls taken from IMT.EVAL Table 1')
end

switch Scenario

    case {'A2'} % InH
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 3*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        if NumLoSConnectionLinks
            loss(LoSConnectionLinks) = iterpar.LoS.PL_A*log10(MsBsDistance(LoSConnectionLinks)) + iterpar.LoS.PL_B + iterpar.LoS.PL_C*log10(CenterFrequency);
            %SF_sigma(LoSConnectionLinks) = 3;
        end

        if NumNLoSConnectionLinks
            loss(NLoSConnectionLinks) = iterpar.NLoS.PL_A*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B + iterpar.NLoS.PL_C*log10(CenterFrequency);
            %SF_sigma(NLoSConnectionLinks) = 4;
        end



    case {'B1'} %UMi        % scenario B1 with d1 & d2 PL model


        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))

            if linkpar.LayoutType % Manhattan
                BsHeight = 10*ones(1,NumLinks);
                MsHeight = 1.5*ones(1,NumLinks);
            else % hexagonal layout
                BsHeight = 10*ones(1,NumLinks);
                MsHeight = 1.5*ones(1,NumLinks);
            end
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        if NumLoSConnectionLinks
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            Dist1 = MsBsDistance_LoS;
            H_bs = BsHeight(LoSConnectionLinks)-1; %effective environment height
            H_ms = MsHeight(LoSConnectionLinks)-1; %effective environment height
            PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
            %SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (Dist1 <= PL_bp & Dist1 >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(Dist1(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C*log10(CenterFrequency);
            end
            %SF_sigma_tmp(ind1) = 3;

            ind2 = (Dist1 >= PL_bp & Dist1 <= iterpar.LoS.PL_range(2));
            if sum(ind2)>0
                loss_LoS(ind2) = iterpar.LoS.PL_A(2)*log10(Dist1(ind2)) + iterpar.LoS.PL_B(2) - 18*log10(H_bs(ind2)) - 18*log10(H_ms(ind2)) + 2*log10(CenterFrequency);
            end
            %SF_sigma_tmp(ind2) = 3;

            %SF_sigma(LoSConnectionLinks) = SF_sigma_tmp;
            loss(LoSConnectionLinks) = loss_LoS;

        end

        if NumNLoSConnectionLinks

            if linkpar.LayoutType % Manhattan

                MsBsDistance_NLoS = MsBsDistance(NLoSConnectionLinks);
                H_bs = BsHeight(NLoSConnectionLinks)-1; %effective environment height
                H_ms = MsHeight(NLoSConnectionLinks)-1; %effective environment height
                PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
                StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
                StreetWidth = StreetWidth(NLoSConnectionLinks);
                loss_LoS = [];

                if isnan(linkpar.Dist1) % NaN default -> distances will be drawn randomly
                    Dist1 = 1;
                    while any(Dist1>MsBsDistance_NLoS-StreetWidth/2 | Dist1<StreetWidth/2)
                        Dist1 = MsBsDistance_NLoS-StreetWidth/2 - (MsBsDistance_NLoS-StreetWidth).*rand(1,length(MsBsDistance_NLoS));
                    end
                else %Dist1 is defined by the user
                    Dist1 = linkpar.Dist1(iterpar.UserIndeces);
                    Dist1 = Dist1(NLoSConnectionLinks);
                    for linkNum = 1:length(MsBsDistance_NLoS) % check applicability and change Dist1 if needed
                        if Dist1(linkNum) > (MsBsDistance_NLoS(linkNum)-StreetWidth(linkNum)/2) || Dist1(linkNum) < StreetWidth(linkNum)/2
                            Dist1(linkNum) = MsBsDistance_NLoS-StreetWidth/2 - (MsBsDistance_NLoS-StreetWidth).*rand(1);
                        end
                    end
                end
                Dist2 = MsBsDistance_NLoS - Dist1;

                loss_a = B1_NLOS_PL(Dist1,Dist2,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);
                loss_b = B1_NLOS_PL(Dist2,Dist1,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);
                loss(NLoSConnectionLinks) = min(loss_a(NLoSConnectionLinks),loss_b(NLoSConnectionLinks));

            else % Hexagonal layout

                loss(NLoSConnectionLinks) = iterpar.NLoS.PL_A_hex*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B_hex+ iterpar.LoS.PL_C_hex*log10(CenterFrequency);

            end

            %SF_sigma(NLoSConnectionLinks) = 4;
        end


    case {'C1'} %SMa
        %Check antenna heights

        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 35*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
        PL_bp = 2*pi*BsHeight.*MsHeight*wimpar.CenterFrequency/2.998e8;

        if NumLoSConnectionLinks %LOS
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            BsHeight_LoS = BsHeight(LoSConnectionLinks);
            MsHeight_LoS = MsHeight(LoSConnectionLinks);
            if isnan(linkpar.BuildingHeight(iterpar.UserIndeces))
                BuildingHeight = repmat(iterpar.LoS.BuildingHeight,1,length(iterpar.UserIndeces));
            else
                BuildingHeight = linkpar.BuildingHeight(iterpar.UserIndeces);
            end
            %SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = 20*log10(40*pi*MsBsDistance_LoS(ind1)*CenterFrequency/3)+min(0.03*BuildingHeight(LoSConnectionLinks).^1.72,10).*log10(MsBsDistance_LoS(ind1))-min(0.044*BuildingHeight(LoSConnectionLinks).^1.72,14.77)+0.002*log10(BuildingHeight(LoSConnectionLinks)).*MsBsDistance_LoS(ind1);
            end
            %SF_sigma_tmp(ind1) = 4;
            ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
            if sum(ind2)>0
                loss_LoS_bp = 20*log10(40*pi*PL_bp*CenterFrequency/3)+min(0.03*BuildingHeight(LoSConnectionLinks).^1.72,10)*log10(PL_bp)-min(0.044*BuildingHeight(LoSConnectionLinks).^1.72,14.77)+0.002*log10(BuildingHeight(LoSConnectionLinks))*PL_bp;
                loss_LoS(ind2) = loss_LoS_bp + 40*log10(MsBsDistance_LoS(ind2)/PL_bp);
            end
            %SF_sigma_tmp(ind2) = 6;
            loss(LoSConnectionLinks) = loss_LoS;
            loss(LoS02ILinks) = loss(LoS02ILinks) + 20;
            loss(LoS02VLinks) = loss(LoS02VLinks) + 9;
            %SF_sigma1(LoSConnectionLinks) = SF_sigma_tmp;
            %SF_sigma2(LosConnectionLinks) = 5; % additional shadowing due to vehicle
        end

        if NumNLoSConnectionLinks %NLOS
            if isnan(linkpar.BuildingHeight(iterpar.UserIndeces))
                BuildingHeight = repmat(iterpar.NLoS.BuildingHeight,1,NumNLoSConnectionLinks);
            else
                BuildingHeight = linkpar.BuildingHeight(iterpar.UserIndeces);
            end
            loss(NLoSConnectionLinks) = 161.04-7.1*log10(StreetWidth(NLoSConnectionLinks))+7.5*log10(BuildingHeight)-(24.37-3.7*(BuildingHeight./BsHeight(NLoSConnectionLinks)).^2).*log10(BsHeight(NLoSConnectionLinks))+(43.42-3.1*log10(BsHeight(NLoSConnectionLinks))).*(log10(MsBsDistance(NLoSConnectionLinks))-3)+20*log10(CenterFrequency)-(3.2*(log10(11.75*MsHeight(NLoSConnectionLinks))).^2-4.97);
            loss(NLoSO2ILinks) = loss(NLoSO2ILinks) + 20; % 20 is wall penetration loss
            loss(NLoSO2VLinks) = loss(NLoSO2VLinks) + 9; % 9 is vehicle penetration loss
            %SF_sigma1(NLoSConnectionLinks) = 8;
            %SF_sigma2(NLosConnectionLinks) = 5; % additional shadowing due to vehicle
        end


    case {'C2'} % UMa
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 25*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
        PL_bp = 4*(BsHeight-1).*(MsHeight-1)*wimpar.CenterFrequency/2.998e8;

        if NumLoSConnectionLinks %LOS
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            BsHeight_LoS = BsHeight(LoSConnectionLinks)-1; %effective environment height
            MsHeight_LoS = MsHeight(LoSConnectionLinks)-1; %effective environment height
            %SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(MsBsDistance_LoS(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C(1)*log10(CenterFrequency);
            end
            %SF_sigma_tmp(ind1) = 4;
            ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
            if sum(ind2)>0
                loss_LoS(ind2) = iterpar.LoS.PL_A(2)*log10(MsBsDistance_LoS(ind2)) + iterpar.LoS.PL_B(2) - 18*log10(MsHeight_LoS(ind2)) - 18*log10(BsHeight_LoS(ind2)) + iterpar.LoS.PL_C(2)*log10(CenterFrequency);
            end
            %SF_sigma_tmp(ind2) = 4;

            loss(LoSConnectionLinks) = loss_LoS + 9; % 9 is vehicle penetration loss
            %SF_sigma1(LoSConnectionLinks) = SF_sigma_tmp;
            %SF_sigma2(LosConnectionLinks) = 5; % additional shadowing due to vehicle

        end

        if NumNLoSConnectionLinks %NLOS
            if isnan(linkpar.BuildingHeight(iterpar.UserIndeces))
                BuildingHeight = repmat(iterpar.NLoS.BuildingHeight,1,NumNLoSConnectionLinks);
            else
                BuildingHeight = linkpar.BuildingHeight(iterpar.UserIndeces);
            end
            loss(NLoSConnectionLinks) = 161.04-7.1*log10(StreetWidth(NLoSConnectionLinks))+7.5*log10(BuildingHeight)-(24.37-3.7*(BuildingHeight./BsHeight(NLoSConnectionLinks)).^2).*log10(BsHeight(NLoSConnectionLinks))+(43.42-3.1*log10(BsHeight(NLoSConnectionLinks))).*(log10(MsBsDistance(NLoSConnectionLinks))-3)+20*log10(CenterFrequency)-(3.2*(log10(11.75*MsHeight(NLoSConnectionLinks))).^2-4.97);
            loss(NLoSConnectionLinks) = loss(NLoSConnectionLinks) + 9; % 9 is vehicle penetration loss
            %SF_sigma1(NLoSConnectionLinks) = 6;
            %SF_sigma2(NLosConnectionLinks) = 5; % additional shadowing due to vehicle
        end



    case {'B4'} % = B1 outdoor-to-indoor
        % Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 10*ones(1,NumLinks);
            MsHeight = (1.5 + 3*(NumFloors-1)).*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
        
        %{
            if isnan(linkpar.Dist1) % NaN default -> will be drawn randomly
            Dist_out = 1;
            while any(Dist_out > MsBsDistance | Dist_out < StreetWidth/2)
                Dist_out = MsBsDistance - (MsBsDistance-StreetWidth/2).*rand(1,length(MsBsDistance));
            end

        else %Dist1 is defined by the user
            Dist_out = linkpar.Dist1(iterpar.UserIndeces);
            Dist_out = Dist_out(NLoSConnectionLinks);
            for linkNum = 1:length(MsBsDistance) % check applicability and change Dist1 if needed
                if Dist_out(linkNum) > MsBsDistance(linkNum) || Dist_out(linkNum) < StreetWidth(linkNum)/2
                    Dist_out(linkNum) = MsBsDistance(linkNum) - (MsBsDistance(linkNum)-StreetWidth(linkNum)/2).*rand(1);
                end
            end
        end
        %}
        tmp = find(MsBsDistance > 25);
        tmp2 = find(MsBsDistance <= 25);
        Dist_in(tmp) = 25*rand(1,length(MsBsDistance(tmp)));  %indoor distance
        Dist_in(tmp2) = MsBsDistance(tmp2).*rand(1,length(MsBsDistance(tmp2)));  %indoor distance
        Dist_out = MsBsDistance - Dist_in;
        Theta = acos(StreetWidth./2./Dist_out)*180/pi; %angle from the BS to the normal of the wall


        %Outdoor loss
        if linkpar.OtoI_OutdoorPL % LOS propagation condition for Outdoor-to-Indoor scenario
            
            % Note!  % LoS --> NLoS since OtoI_OutdoorPL is a subscenario
            MsBsDistance_LoS = MsBsDistance(NLoSConnectionLinks); 
            Dist1 = MsBsDistance_LoS;
            H_bs = BsHeight(NLoSConnectionLinks)-1; %effective environment height
            H_ms = MsHeight(NLoSConnectionLinks)-1; %effective environment height
            PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
            %SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (Dist1 <= PL_bp & Dist1 >= fixpar.B1.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = fixpar.B1.LoS.PL_A(1)*log10(Dist1(ind1)) + fixpar.B1.LoS.PL_B(1) + fixpar.B1.LoS.PL_C*log10(CenterFrequency);
            end
            %SF_sigma_tmp(ind1) = 3;

            ind2 = (Dist1 >= PL_bp & Dist1 <= fixpar.B1.LoS.PL_range(2));
            if sum(ind2)>0
                loss_LoS(ind2) = fixpar.B1.LoS.PL_A(2)*log10(Dist1(ind2)) + fixpar.B1.LoS.PL_B(2) - 18*log10(H_bs(ind2)) - 18*log10(H_ms(ind2)) + 2*log10(CenterFrequency);
            end
            %SF_sigma_tmp(ind2) = 3;

            %SF_sigma(NLoSConnectionLinks) = SF_sigma_tmp;
            loss(NLoSConnectionLinks) = loss_LoS;

        else % NLOS

            if linkpar.LayoutType % Manhattan

                MsBsDistance_NLoS = MsBsDistance(NLoSConnectionLinks);
                H_bs = BsHeight(NLoSConnectionLinks)-1; %effective environment height
                H_ms = MsHeight(NLoSConnectionLinks)-1; %effective environment height
                PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
                StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
                StreetWidth = StreetWidth(NLoSConnectionLinks);
                loss_LoS = [];

                if isnan(linkpar.Dist1) % NaN default -> distances will be drawn randomly
                    Dist1 = 1;
                    while any(Dist1>MsBsDistance_NLoS-StreetWidth/2 | Dist1<StreetWidth/2)
                        Dist1 = MsBsDistance_NLoS-StreetWidth/2 - (MsBsDistance_NLoS-StreetWidth).*rand(1,length(MsBsDistance_NLoS));
                    end
                else %Dist1 is defined by the user
                    Dist1 = linkpar.Dist1(iterpar.UserIndeces);
                    Dist1 = Dist1(NLoSConnectionLinks);
                    for linkNum = 1:length(MsBsDistance_NLoS) % check applicability and change Dist1 if needed
                        if Dist1(linkNum) > (MsBsDistance_NLoS(linkNum)-StreetWidth(linkNum)/2) || Dist1(linkNum) < StreetWidth(linkNum)/2
                            Dist1(linkNum) = MsBsDistance_NLoS-StreetWidth/2 - (MsBsDistance_NLoS-StreetWidth).*rand(1);
                        end
                    end
                end
                Dist2 = MsBsDistance_NLoS - Dist1;

                loss_a = B1_NLOS_PL(Dist1,Dist2,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);
                loss_b = B1_NLOS_PL(Dist2,Dist1,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);
                loss(NLoSConnectionLinks) = min(loss_a(NLoSConnectionLinks),loss_b(NLoSConnectionLinks));

            else % Hexagonal layout

                % loss(NLoSConnectionLinks) = iterpar.NLoS.PL_A_hex*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B_hex+ iterpar.LoS.PL_C_hex*log10(CenterFrequency);
                loss(NLoSConnectionLinks) = 36.7*log10(MsBsDistance(NLoSConnectionLinks)) + 22.7+ 26*log10(CenterFrequency);
            end

            %SF_sigma(NLoSConnectionLinks) = 4;
        end

        LossOut = loss;
        %indoor loss
        LossIn = 0.5*Dist_in; %alpha is 0.5 dB/meter
        %Through wall loss
        if linkpar.LayoutType % Manhattan
            LossWall = 14+15*(1-cos(Theta)).^2;
        else % Hexagonal
            LossWall = 20;
        end
        %Total loss
        loss = LossOut+LossIn+LossWall;

        %SF_sigma(NLoSConnectionLinks) = 7;



    case {'D1'} % RMa
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 35*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
        PL_bp = 2*pi*BsHeight.*MsHeight*wimpar.CenterFrequency/2.998e8;
        
        if NumLoSConnectionLinks %LOS
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            BsHeight_LoS = BsHeight(LoSConnectionLinks);
            MsHeight_LoS = MsHeight(LoSConnectionLinks);
            if isnan(linkpar.BuildingHeight(iterpar.UserIndeces))
                BuildingHeight = repmat(iterpar.LoS.BuildingHeight,1,NumLoSConnectionLinks);
            else
                BuildingHeight = linkpar.BuildingHeight(iterpar.UserIndeces);
            end
            %SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = 20*log10(40*pi*MsBsDistance_LoS(ind1)*CenterFrequency/3)+min(0.03*BuildingHeight.^1.72,10).*log10(MsBsDistance_LoS(ind1))-min(0.044*BuildingHeight.^1.72,14.77)+0.002*log10(BuildingHeight).*MsBsDistance_LoS(ind1);
            end
            %SF_sigma_tmp(ind1) = 4;
            ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
            if sum(ind2)>0
                %loss_LoS_bp = 20*log10(40*pi*PL_bp*CenterFrequency/3)+min(0.03*h^1.72,10)*log10(PL_bp)-min(0.044*h^1.72,14.77)+0.002*log10(h)*PL_bp;
                loss_LoS_bp = 20*log10(40*pi*PL_bp*CenterFrequency/3)+min(0.03*BuildingHeight^1.72,10)*log10(PL_bp)-min(0.044*BuildingHeight^1.72,14.77)+0.002*log10(BuildingHeight)*PL_bp;
                loss_LoS(ind2) = loss_LoS_bp + 40*log10(MsBsDistance_LoS(ind2)/PL_bp);
            end
            %SF_sigma_tmp(ind2) = 6;
            
            loss(LoSConnectionLinks) = loss_LoS + 9; % 9 is vehicle penetration loss
            %SF_sigma1(LoSConnectionLinks) = SF_sigma_tmp;
            %SF_sigma2(LosConnectionLinks) = 5; % additional shadowing due to vehicle
        end

        if NumNLoSConnectionLinks %NLOS
            if isnan(linkpar.BuildingHeight(iterpar.UserIndeces))
                BuildingHeight = repmat(iterpar.NLoS.BuildingHeight,1,NumLinks);
            else
                BuildingHeight = linkpar.BuildingHeight(iterpar.UserIndeces);
            end
            loss(NLoSConnectionLinks) = 161.04-7.1*log10(StreetWidth(NLoSConnectionLinks))+7.5*log10(BuildingHeight(NLoSConnectionLinks))-(24.37-3.7*(BuildingHeight(NLoSConnectionLinks)./BsHeight(NLoSConnectionLinks)).^2).*log10(BsHeight(NLoSConnectionLinks))+(43.42-3.1*log10(BsHeight(NLoSConnectionLinks))).*(log10(MsBsDistance(NLoSConnectionLinks))-3)+20*log10(CenterFrequency)-(3.2*(log10(11.75*MsHeight(NLoSConnectionLinks))).^2-4.97);
            loss(NLoSConnectionLinks) = loss(NLoSConnectionLinks) + 9; % 9 is vehicle penetration loss
            %SF_sigma1(NLoSConnectionLinks) = 6;
            %SF_sigma2(NLosConnectionLinks) = 5; % additional shadowing due to vehicle
        end


    otherwise       % all other scenarios with one d PL model

        BsHeight = 10*ones(1,NumLinks);
        MsHeight = 1.5*ones(1,NumLinks);
        %SF_sigma = 3*ones(1,NumLinks);
        loss(LoSConnectionLinks) = iterpar.LoS.PL_B + iterpar.LoS.PL_A*log10(MsBsDistance(LoSConnectionLinks));
        loss(NLoSConnectionLinks) = iterpar.NLoS.PL_B + iterpar.NLoS.PL_A*log10(MsBsDistance(NLoSConnectionLinks));

end     % end switch Scenario


% output
linkpar.MsHeight(iterpar.UserIndeces) = MsHeight;
linkpar.BsHeight(iterpar.UserIndeces) = BsHeight;
if ~isempty(LoS02ILinks); linkpar.LoS02ILinks(iterpar.UserIndeces(LoS02ILinks)) = LoS02ILinks; end
if ~isempty(LoS02VLinks); linkpar.LoS02VLinks(iterpar.UserIndeces(LoS02VLinks)) = LoS02VLinks; end
if ~isempty(NLoSO2ILinks); linkpar.NLoS02ILinks(iterpar.UserIndeces(NLoSO2ILinks)) = NLoSO2ILinks; end
if ~isempty(NLoSO2VLinks); linkpar.NLoS02VLinks(iterpar.UserIndeces(NLoSO2VLinks)) = NLoSO2VLinks; end
%iterpar.NLoS.SF_sigma=SF_sigma(1);
loss=loss(:);

% function d=distrnd(num,rmax)
% % DISTRND Distance from BS in a circular cell
% %   D=DISTRND(K,RMAX) generates K random variables from the pdf
% %   p(r)=2*r/RMAX^2. This is the pdf of distance from base station when
% %   users are uniformly (in area) distributed in a cell with radius RMAX
% %   1 x num vector.
%
% %   Authors: Jari Salo (HUT), Marko Milojevic (TUI)
%
% % create random variables from triangular pdf whose width is 2*rmax
% a=sum(repmat(rmax,2,1).*rand(2,num));
%
% % fold the random variables about the rmax
% inds=find(a>rmax);
% a(inds)=-a(inds)+2*rmax(inds);
%
% d=a(:).';


function loss = B1_NLOS_PL(Dist1,Dist2,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);


% Same as B1 LOS
ind1 = (Dist1 <= PL_bp & Dist1 >= (iterpar.NLoS.PL_range(1)-StreetWidth/2));
if sum(ind1)>0
    loss_LoS(ind1) = iterpar.NLoS.PL_A(1)*log10(Dist1(ind1)) + iterpar.NLoS.PL_B(1) + 20*log10(CenterFrequency);
end

ind2 = (Dist1 >= PL_bp & Dist1 <= iterpar.NLoS.PL_range(2)-StreetWidth/2);
if sum(ind2)>0
    loss_LoS(ind2) = iterpar.NLoS.PL_A(2)*log10(Dist1(ind2)) + iterpar.NLoS.PL_B(2) - 18*log10(H_bs(ind2)) - 18*log10(H_ms(ind2)) + 2*log10(CenterFrequency);
end

% plus component for B1 NLOS
nj = max(2.8-0.0024*Dist1 , 1.84);
loss(NLoSConnectionLinks) = loss_LoS + 17.9 - 12.5*nj + 10*nj.*log10(Dist2) + 3*log10(CenterFrequency);
