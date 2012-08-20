%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed delays,  powers and cluster-wice K-factors
%  for different scenarios.
% Needed when wimpar.FixedPdpUsed='yes'
function [taus_LoS,Pprime_LoS,Kcluster_LoS, ...
        taus_NLoS,Pprime_NLoS,Kcluster_NLoS] = fixedPdp(wimpar,iterpar)

LoSConnectionLinks = iterpar.LoSConnectionLinks;
NLoSConnectionLinks = iterpar.NLoSConnectionLinks;
Scenario = iterpar.Scenario;


switch Scenario
   
    case {'A2'}     % InH
        taus_LoS = [0 10 25 25 30 35 40 50 55 60 65 85 90 100 130]*1E-9;
        Pprime_LoS = 10.^(-[0 15.7 10.5 16.7 17.6 14.1 12.9 19.5 21.8 20.8 24.1 13.9 20.1 18.0 21.0]/10);
        Kcluster_LoS = [5 ;...               % K-factors for CDL clusters [dB]
                       1];                   % cluster number
                   
        taus_NLoS = [0 15 20 25 30 40 55 60 60 70 75 75 90 150 160 170 195 205 225]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[2.4 1.9 8.1 1.8 0 2.3 3.7 8.4 3.2 9.7 6.2 8.9 4.0 14.1 12.1 10.6 19.6 16.8 13.5]/10);          % lin.
        Kcluster_NLoS = [-10000000  ;...     % K-factors for CDL clusters [dB]
                       1   ];                % cluster number
        
    case {'B1'}     % UMi
        taus_LoS = [0 30 85 135 160 195 210 255 280 340 360 420]*1E-9;
        Pprime_LoS = 10.^(-[0 15.6 14.0 15.8 20.7 17.3 21.8 17.7 21.6 23.0 24.5 25.0]/10);
        Kcluster_LoS = [6.0 ;...             % K-factors for CDL clusters [dB]
                       1];                   % cluster number

        %this one is not updated accordingly in D111, 15.2.2007. Values from LH
        taus_NLoS = [0 10 20 35 40 55 55 200 205 250 330 440 440 515 530 580 590 625 730]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[6.7 4.9 1.9 6.3 3 7.5 6.4 10.8 5.2 4.9 9.2 15.5 16.7 12.4 16.9 12.7 23.5 22.1 23.6]/10); % lin.
        Kcluster_NLoS = [-100000 ;...        % K-factors for CDL clusters [dB]
                      1 ];                   % cluster number
     
      case {'B4'}   % UMi O-to-I 

        taus_NLoS = [0 0 5 10 35 35 65 120 125 195 250 305]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[13.0 21.7 16.7 24.9 29.2 19.9 13.4 23.3 33.7 29.1 34.0 35.9]/10);          % lin.
        Kcluster_NLoS = [-10000000  ;...            % K-factors for CDL clusters [dB]
                       1   ];                       % cluster number
                 
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;
        
    case {'C1'}     % SMa

        taus_LoS = [0 85 135 135 170 190 275 290 290 410 445 500 620 655 960]*1E-9;
        Pprime_LoS = 10.^(-[0 21.6 26.3 25.1 25.4 22.0 29.2 24.3 23.2 32.2 26.5 32.1 28.5 30.5 32.6]/10);
        Kcluster_LoS = [12.9;...                   % K-factors for CDL clusters [dB]
                     1 ];                      % cluster number

        taus_NLoS = [0 25 35 35 45 65 65 75 145 160 195 200 205 770]*1E-9;                                % [s]
        Pprime_NLoS = 10.^(-[3.0 7.5 10.5 3.2 6.1 14.0 6.4 3.1 4.6 8.0 7.2 3.1 9.5 22.4]/10); % lin.
        Kcluster_NLoS = [-10000000;...                    % K-factors for CDL clusters [dB]
                            1   ];                   % cluster number

    case {'C2'}     % UMa
        
        taus_LoS = [0 15 30 45 220 310 365 440 450 535 595 640]*1E-9;
        Pprime_LoS = 10.^(-[0 15.4 12.6 14.1 19.4 23.8 16.7 19.4 24.8 21.3 23.9 25.0]/10);
        Kcluster_LoS = [7.0;...                   % K-factors for CDL clusters [dB]
                     1 ];                      % cluster number
        
        taus_NLoS = [0 5 20 45 265 290 325 340 355 440 555 645 970 1015 1220 1395 1540 1750 1870 1885]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[3.5 9.2 3.0 7.8 3.7 8.6 2.5 7.3 3.8 6.9 8.9 9.0 9.8 15.0 13.4 14.9 16.7 11.2 18.2 17.8]/10);                      % lin.
        Kcluster_NLoS = [-1000000;...                          % K-factors for CDL clusters [dB]
                           1];                            % cluster number
        

    case {'D1'}    % RMa
        taus_LoS = [0 35 45 65 65 110 125 125 170 170 200]*1E-9;
        Pprime_LoS = 10.^(-[16.8 18.3 21.2 17.1 19.7 23.8 22.9 20.9 21.9]/10);
        Kcluster_LoS = [13.7 ;...                   % K-factors for CDL clusters [dB]
                       1 ];                      % cluster number

        taus_NLoS = [0 0 5 10 20 15 55 100 110 220]*1E-9;                                % [s]
        Pprime_NLoS = 10.^(-[4.9 7.8 5.2 2.7 2.6 5.8 2.7 5.6 7.3 10.3]/10); % lin.
        Kcluster_NLoS = [-10000000;...                    % K-factors for CDL clusters [dB]
                      1   ];                      % cluster number

end % switch