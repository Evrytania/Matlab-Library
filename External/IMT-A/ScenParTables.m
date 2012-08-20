
function fixpar = ScenParTables(StreetWidth)
%   SCENPARTABLES Scenario specific parameters of WIM channel model
%   FIXPAR = SCENPARTABLES(StreetWidth) Sets WIM scenario
%   specific parameters given in [1, Table  A1-7 and A1-2]. These
%   parameters are used in GENERATE_BULK_PAR to generate spatio-temporal
%   channel parameters for channel coefficient generation.
%
%   Ref. [1]: ITU-R IMT.EVAL Channel model, Finland contribution Sept 2008
%
%   See also WIMPARSET and GENERATE_BULK_PAR.

%   Authors: Pekka Ky�sti (EBIT), Lassi Hentil� (EBIT), Daniela Laselva
%   (EBIT), Marko Milojevic (TUI), Mikko Alatossava (CWC/UOULU),Jianhua
%   Zhang(BUPT),Guangyi Liu(CMCC)
%

%   Modifications:
%   PL_bp parameter removed. Now calculated in pathloss.    15.5.2006 PekKy
%   A1 lambda parameters updated.                           15.5.2006 PekKy
%   B1 path loss parameters updated.                        8.9.2006  PekKy
%   update of proportionality factors                       5.10.2006 HentLas
%   D1.1.1 parameters - old updated, new introduced         15.1.2007 Marko
%   Output includes all scenarios and LoS and NLoS          12.2.2007 MikkoA
%   Angle proportionality factors removed,... 
%   K-factor definition changed,...
%   Path loss frequency dependence added,...
%   Scenario C4 added,...
%   Cross-corr. and decorr.dist added for K-factor          3.10.2007 HentLas
%   Parameters updated to comply with ITU-R IMT.EVAL ver    
%    Note! Path loss parameters still ambiguous, check      26.5.2008 PekKy
%   Parameters updated to match with Seoul contribution     19.9.2008 HentLas, PekKy


%Dispersion parameters [1, Table 2]

%% InH A2, LoS
% Fixed scenario specific parameters
fixpar.A2.LoS.NumClusters = 15;         % Number of ZDSC    [1, Table 2]
fixpar.A2.LoS.r_DS   = 3.6;             % delays spread proportionality factor
fixpar.A2.LoS.PerClusterAS_D = 5;       % Per cluster FS angle spread [deg] [1, Table 2]
fixpar.A2.LoS.PerClusterAS_A = 8;       % Per cluster MS angle spread [deg] [1, Table 2]
fixpar.A2.LoS.LNS_ksi = 6;              % ZDSC LNS ksi [dB], per cluster shadowing [1, Table 2]

% Cross correlation coefficients [1, Table 2]
fixpar.A2.LoS.asD_ds = 0.6;             % departure AS vs delay spread
fixpar.A2.LoS.asA_ds = 0.8;             % arrival AS vs delay spread
fixpar.A2.LoS.asA_sf = -0.5;            % arrival AS vs shadowing std
fixpar.A2.LoS.asD_sf = -0.4;            % departure AS vs shadowing std
fixpar.A2.LoS.ds_sf  = -0.8;            % delay spread vs shadowing std
fixpar.A2.LoS.asD_asA = 0.4;            % departure AS vs arrival AS
fixpar.A2.LoS.asD_kf = 0;               % departure AS vs k-factor
fixpar.A2.LoS.asA_kf = 0;               % arrival AS vs k-factor
fixpar.A2.LoS.ds_kf = -0.5;             % delay spread vs k-factor
fixpar.A2.LoS.sf_kf = 0.5;              % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.A2.LoS.xpr_mu    = 11;           % XPR mean [dB]
fixpar.A2.LoS.xpr_sigma = 0;            % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.A2.LoS.DS_mu      = -7.70;       % delay spread, mean [log10(s)]
fixpar.A2.LoS.DS_sigma   = 0.18;        % delay spread, std [log10(s)]
fixpar.A2.LoS.AS_D_mu    = 1.60;        % arrival angle spread, mean [log10(deg)]
fixpar.A2.LoS.AS_D_sigma = 0.18;        % arrival angle spread, std [log10(deg)]
fixpar.A2.LoS.AS_A_mu    = 1.62;        % departure angle spread, mean [log10(deg)]
fixpar.A2.LoS.AS_A_sigma = 0.22;        % departure angle spread, std [log10(deg)]
fixpar.A2.LoS.SF_sigma   = 3;           % shadowing std [dB] (zero mean)
fixpar.A2.LoS.KF_mu      = 7;           % K-factor mean [dB]
fixpar.A2.LoS.KF_sigma   = 4;           % K-factor std [dB]

% "Decorrelation distances" [1, Table 2]
fixpar.A2.LoS.DS_lambda   = 8;          % [m], delay spread
fixpar.A2.LoS.AS_D_lambda = 7;          % [m], departure azimuth spread
fixpar.A2.LoS.AS_A_lambda = 5;          % [m], arrival azimuth spread
fixpar.A2.LoS.SF_lambda   = 10;         % [m], shadowing
fixpar.A2.LoS.KF_lambda   = 4;          % [m], k-factor  

% Path loss PL = Alog10(d) + B + Clog10(fc)  [1, Table 1]
fixpar.A2.LoS.PL_A = 16.9;              % path loss exponent
fixpar.A2.LoS.PL_B = 32.8;              % path loss intercept
fixpar.A2.LoS.PL_C = 20;                % path loss frequency dependence factor
fixpar.A2.LoS.PL_range = [3 100];       % applicability range [m], (min max)
%%


%% InH A2, NLoS
% Fixed scenario specific parameters
fixpar.A2.NLoS.NumClusters = 19;        % Number of ZDSC    [1, Table 2]
fixpar.A2.NLoS.r_DS   = 3.0;            % delays spread proportionality factor
fixpar.A2.NLoS.PerClusterAS_D = 5;      % Per cluster FS angle spread [deg] [1, Table 2]
fixpar.A2.NLoS.PerClusterAS_A = 11;     % Per cluster MS angle spread [deg] [1, Table 2]
fixpar.A2.NLoS.LNS_ksi = 3;             % ZDSC LNS ksi [dB], per cluster shadowing [1, Table 2]

% Cross correlation coefficients [1, Table 2]
fixpar.A2.NLoS.asD_ds =  0.4;           % departure AS vs delay spread
fixpar.A2.NLoS.asA_ds =  0;             % arrival AS vs delay spread
fixpar.A2.NLoS.asA_sf =  -0.4;          % arrival AS vs shadowing std
fixpar.A2.NLoS.asD_sf =  0;             % departure AS vs shadowing std
fixpar.A2.NLoS.ds_sf  =  -0.5;          % delay spread vs shadowing std
fixpar.A2.NLoS.asD_asA = 0;             % departure AS vs arrival AS
fixpar.A2.NLoS.asD_kf = 0;              % departure AS vs k-factor
fixpar.A2.NLoS.asA_kf = 0;              % arrival AS vs k-factor
fixpar.A2.NLoS.ds_kf  = 0;              % delay spread vs k-factor
fixpar.A2.NLoS.sf_kf  = 0;              % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.A2.NLoS.xpr_mu    = 10;          % XPR mean [dB]
fixpar.A2.NLoS.xpr_sigma = 0;           % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.A2.NLoS.DS_mu      = -7.41;      % delay spread, mean [log10(s)]
fixpar.A2.NLoS.DS_sigma   = 0.14;       % delay spread, std [log10(s)]
fixpar.A2.NLoS.AS_D_mu    = 1.62;       % arrival angle spread, mean [log10(deg)]
fixpar.A2.NLoS.AS_D_sigma = 0.25;       % arrival angle spread, std [log10(deg)]
fixpar.A2.NLoS.AS_A_mu    = 1.77;       % departure angle spread, mean [log10(deg)]
fixpar.A2.NLoS.AS_A_sigma = 0.16;       % departure angle spread, std [log10(deg)]
fixpar.A2.NLoS.SF_sigma   = 4;          % shadowing std [dB] (zero mean)
fixpar.A2.NLoS.KF_mu      = 0;          % k-factor, dummy value
fixpar.A2.NLoS.KF_sigma   = 0;          % k-factor, dummy value

% "Decorrelation distances" [1, Table 2]
fixpar.A2.NLoS.DS_lambda   = 5;         % [m], delay spread
fixpar.A2.NLoS.AS_D_lambda = 3;         % [m], departure azimuth spread
fixpar.A2.NLoS.AS_A_lambda = 3;         % [m], arrival azimuth spread
fixpar.A2.NLoS.SF_lambda   = 6;         % [m], shadowing
fixpar.A2.NLoS.KF_lambda   = NaN;       % [m], k-factor  

% Path loss PL = Alog10(d) + B + Clog10(fc)  [1, Table 1]
fixpar.A2.NLoS.PL_A = 43.3;             % path loss exponent
fixpar.A2.NLoS.PL_B = 11.5;             % path loss intercept
fixpar.A2.NLoS.PL_C = 20;               % path loss frequency dependence factor
fixpar.A2.NLoS.PL_X = [NaN NaN];        % path loss wall factor [1, Table 1]
fixpar.A2.NLoS.PL_range = [10 150];     % applicability range [m], (min max)
%%


%% UMi B1, LoS
% Fixed scenario specific parameters
fixpar.B1.LoS.NumClusters = 12;         % Number of ZDSC    [1, Table 2]
fixpar.B1.LoS.r_DS   = 3.2;             % delays spread proportionality factor
fixpar.B1.LoS.PerClusterAS_D = 3;       % Per cluster FS angle spread [deg] [1, Table 2]
fixpar.B1.LoS.PerClusterAS_A = 17;      % Per cluster MS angle spread [deg] [1, Table 2]
fixpar.B1.LoS.LNS_ksi = 3;              % ZDSC LNS ksi [dB], per cluster shadowing [1, Table 2]

% Cross correlation coefficients [1, Table 2]
fixpar.B1.LoS.asD_ds = 0.5;             % departure AS vs delay spread
fixpar.B1.LoS.asA_ds = 0.8;             % arrival AS vs delay spread
fixpar.B1.LoS.asA_sf = -0.4;            % arrival AS vs shadowing std
fixpar.B1.LoS.asD_sf = -0.5;            % departure AS vs shadowing std
fixpar.B1.LoS.ds_sf  = -0.4;            % delay spread vs shadowing std
fixpar.B1.LoS.asD_asA = 0.4;            % departure AS vs arrival AS
fixpar.B1.LoS.asD_kf = -0.2;            % departure AS vs k-factor, PK 18.8.08
fixpar.B1.LoS.asA_kf = -0.3;            % arrival AS vs k-factor,   PK 18.8.08
fixpar.B1.LoS.ds_kf = -0.7;             % delay spread vs k-factor
fixpar.B1.LoS.sf_kf = 0.5;              % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.B1.LoS.xpr_mu    = 9;            % XPR mean [dB]
fixpar.B1.LoS.xpr_sigma = 0;            % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
fixpar.B1.LoS.DS_mu      = -7.19;       % delay spread, mean [log10(s)]
fixpar.B1.LoS.DS_sigma   = 0.40;        % delay spread, std [log10(s)]
fixpar.B1.LoS.AS_D_mu    = 1.20;        % arrival angle spread, mean [log10(deg)]
fixpar.B1.LoS.AS_D_sigma = 0.43;        % arrival angle spread, std [log10(deg)]
fixpar.B1.LoS.AS_A_mu    = 1.75;        % departure angle spread, mean [log10(deg)]
fixpar.B1.LoS.AS_A_sigma = 0.19;        % departure angle spread, std [log10(deg)]
fixpar.B1.LoS.SF_sigma   = 3;           % shadowing std [dB] (zero mean)
fixpar.B1.LoS.KF_mu = 9;                % K-factor mean [dB]
fixpar.B1.LoS.KF_sigma = 5;             % K-factor std [dB]

% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.B1.LoS.DS_lambda   = 7;          % [m], delay spread
fixpar.B1.LoS.AS_D_lambda = 8;          % [m], departure azimuth spread
fixpar.B1.LoS.AS_A_lambda = 8;          % [m], arrival azimuth spread
fixpar.B1.LoS.SF_lambda   = 10;         % [m], shadowing
fixpar.B1.LoS.KF_lambda   = 15;         % [m], k-factor 

% Path loss PL = Alog10(d) + B + Clog10(fc/5)  [1, Table 1]
fixpar.B1.LoS.PL_A = [22.0 40.0];       % path loss exponent, [d<d_bp d>d_bp]
fixpar.B1.LoS.PL_B = [28.0 7.8];        % path loss intercept, [d<d_bp d>d_bp]
fixpar.B1.LoS.PL_C = 20;                % path loss frequency dependence factor
fixpar.B1.LoS.PL_range = [10 5000];     % applicability range [m], (min max)
%%


%% UMi B1, NLoS
% Fixed scenario specific parameters
fixpar.B1.NLoS.NumClusters = 19;        % Number of ZDSC    [1, Table 2]
fixpar.B1.NLoS.r_DS   = 3;              % delays spread proportionality factor
fixpar.B1.NLoS.PerClusterAS_D = 10;     % Per cluster FS angle spread [deg] [1, Table 2]
fixpar.B1.NLoS.PerClusterAS_A = 22;     % Per cluster MS angle spread [deg] [1, Table 2]
fixpar.B1.NLoS.LNS_ksi = 3;             % ZDSC LNS ksi [dB], per cluster shadowing [1, Table 2]

% Cross correlation coefficients [1, Table 2]
fixpar.B1.NLoS.asD_ds =  0;             % departure AS vs delay spread
fixpar.B1.NLoS.asA_ds =  0.4;           % arrival AS vs delay spread
fixpar.B1.NLoS.asA_sf = -0.4;           % arrival AS vs shadowing std
fixpar.B1.NLoS.asD_sf = 0;              % departure AS vs shadowing std
fixpar.B1.NLoS.ds_sf  = -0.7;           % delay spread vs shadowing std
fixpar.B1.NLoS.asD_asA = 0;             % departure AS vs arrival AS
fixpar.B1.NLoS.asD_kf = 0;              % departure AS vs k-factor
fixpar.B1.NLoS.asA_kf = 0;              % arrival AS vs k-factor
fixpar.B1.NLoS.ds_kf = 0;               % delay spread vs k-factor
fixpar.B1.NLoS.sf_kf = 0;               % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.B1.NLoS.xpr_mu    = 8;           % XPR mean [dB]
fixpar.B1.NLoS.xpr_sigma = 0;           % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.B1.NLoS.DS_mu      = -6.89;      % delay spread, mean [log10(s)]
fixpar.B1.NLoS.DS_sigma   = 0.54;       % delay spread, std [log10(s)]
fixpar.B1.NLoS.AS_D_mu    = 1.41;       % arrival angle spread, mean [log10(deg)]
fixpar.B1.NLoS.AS_D_sigma = 0.17;       % arrival angle spread, std [log10(deg)]
fixpar.B1.NLoS.AS_A_mu    = 1.84;       % departure angle spread, mean [log10(deg)]
fixpar.B1.NLoS.AS_A_sigma = 0.15;       % departure angle spread, std [log10(deg)]
fixpar.B1.NLoS.SF_sigma   = 4;          % shadowing std [dB] (zero mean)
fixpar.B1.NLoS.KF_mu      = 0;          % k-factor, dummy value
fixpar.B1.NLoS.KF_sigma   = 0;          % k-factor, dummy value

%%% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.B1.NLoS.DS_lambda   = 10;        % [m], delay spread
fixpar.B1.NLoS.AS_D_lambda = 10;        % [m], departure azimuth spread
fixpar.B1.NLoS.AS_A_lambda = 9;         % [m], arrival azimuth spread
fixpar.B1.NLoS.SF_lambda   = 13;        % [m], shadowing
fixpar.B1.NLoS.KF_lambda   = NaN;       % [m], k-factor 

% Path loss, Note! see the path loss equation... HEXAGONAL

fixpar.B1.NLoS.PL_A_hex = 36.7;             % path loss exponent
fixpar.B1.NLoS.PL_B_hex = 22.7;             % path loss intercept, 
fixpar.B1.LoS.PL_C_hex = 26;                % path loss frequency dependence factor
fixpar.B1.NLoS.PL_range = [10 2000];        % applicability range [m], (min max)
% B1 PL for MANHATTAN GRID!
fixpar.B1.NLoS.PL_A = [22.0 40.0];          % path loss exponent, [d<d_bp d>d_bp]
fixpar.B1.NLoS.PL_B = [28.0 7.8];           % path loss intercept, [d<d_bp d>d_bp]
fixpar.B1.LoS.PL_C = 20;                    % path loss frequency dependence factor
fixpar.B1.NLoS.PL_range_Manhattan = 5000;   % applicability range upper limit [m]

%%


%% SMa C1, LoS
% Fixed scenario specific parameters [1, Table 2]
fixpar.C1.LoS.NumClusters = 15;     % Number of ZDSC
fixpar.C1.LoS.r_DS   = 2.4;         % delays spread proportionality factor
fixpar.C1.LoS.PerClusterAS_D = 5;   % Per cluster FS angle spread [deg]
fixpar.C1.LoS.PerClusterAS_A = 5;   % Per cluster MS angle spread [deg]
fixpar.C1.LoS.LNS_ksi = 3;          % ZDSC LNS ksi [dB], per cluster shadowing

% Cross correlation coefficients
fixpar.C1.LoS.asD_ds = 0;           % departure AS vs delay spread
fixpar.C1.LoS.asA_ds = 0.8;         % arrival AS vs delay spread
fixpar.C1.LoS.asA_sf = -0.5;        % arrival AS vs shadowing std
fixpar.C1.LoS.asD_sf = -0.5;        % departure AS vs shadowing std
fixpar.C1.LoS.ds_sf  = -0.6;        % delay spread vs shadowing std
fixpar.C1.LoS.asD_asA = 0;          % departure AS vs arrival AS
fixpar.C1.LoS.asD_kf = 0;           % departure AS vs k-factor
fixpar.C1.LoS.asA_kf = 0;           % arrival AS vs k-factor
fixpar.C1.LoS.ds_kf = 0;            % delay spread vs k-factor
fixpar.C1.LoS.sf_kf = 0;            % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.C1.LoS.xpr_mu    = 8;        % XPR mean [dB]
fixpar.C1.LoS.xpr_sigma = 0;        % XPR std  [dB], PK 18.8.2008

%%% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.C1.LoS.DS_mu      = -7.23;   % delay spread, mean [s-dB]
fixpar.C1.LoS.DS_sigma   = 0.38;    % delay spread, std [s-dB]
fixpar.C1.LoS.AS_D_mu    = 0.78;    % arrival angle spread, mean [deg-dB]
fixpar.C1.LoS.AS_D_sigma = 0.12;    % arrival angle spread, std [deg-dB]
fixpar.C1.LoS.AS_A_mu    = 1.48;    % departure angle spread, mean [deg-dB]
fixpar.C1.LoS.AS_A_sigma = 0.20;    % departure angle spread, std [deg-dB]
fixpar.C1.LoS.SF_sigma   = 4;       % shadowing std [dB] (zero mean)
fixpar.C1.LoS.KF_mu = 9;            % K-factor mean [dB]
fixpar.C1.LoS.KF_sigma = 7;         % K-factor std [dB]

%%% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.C1.LoS.DS_lambda   = 6;       % [m], delay spread
fixpar.C1.LoS.AS_D_lambda = 15;      % [m], departure azimuth spread
fixpar.C1.LoS.AS_A_lambda = 20;      % [m], arrival azimuth spread
fixpar.C1.LoS.SF_lambda   = 40;      % [m], shadowing
fixpar.C1.LoS.KF_lambda   = 10;      % [m], k-factor

% Path loss, Note! see the path loss equation...
fixpar.C1.LoS.PL_A = [NaN NaN];         % path loss exponent [dB], [d<d_bp d>d_bp]
fixpar.C1.LoS.PL_B = [NaN NaN];         % path loss intercept [dB], [d<d_bp d>d_bp]
fixpar.C1.LoS.PL_C = [NaN NaN];         % path loss frequency dependence factor [dB], [d<d_bp d>d_bp]
fixpar.C1.LoS.PL_range = [10 5000];     % applicability range [m]
fixpar.C1.LoS.BuildingHeight = 10;      % average building height
%%


%% SMa C1, NLoS
%%% Fixed scenario specific parameters
fixpar.C1.NLoS.NumClusters = 14;       % Number of ZDSC    [1, Table 2]
fixpar.C1.NLoS.r_DS   = 1.5;           % delays spread proportionality factor [1, Table 2]
fixpar.C1.NLoS.PerClusterAS_D = 2;     % Per cluster FS angle spread [deg] [1, Table 2]
fixpar.C1.NLoS.PerClusterAS_A = 10;    % Per cluster MS angle spread [deg] [1, Table 2]
fixpar.C1.NLoS.LNS_ksi = 3;            % ZDSC LNS ksi [dB], per cluster shadowing [1, Table 2]

% Cross correlation coefficients [1, Table 2]
fixpar.C1.NLoS.asD_ds = 0;              % departure AS vs delay spread
fixpar.C1.NLoS.asA_ds = 0.7;            % arrival AS vs delay spread
fixpar.C1.NLoS.asA_sf = 0;              % arrival AS vs shadowing std
fixpar.C1.NLoS.asD_sf = -0.4;           % departure AS vs shadowing std
fixpar.C1.NLoS.ds_sf  = -0.4;           % delay spread vs shadowing std
fixpar.C1.NLoS.asD_asA = 0;             % departure AS vs arrival AS
fixpar.C1.NLoS.asD_kf = 0;              % departure AS vs k-factor
fixpar.C1.NLoS.asA_kf = 0;              % arrival AS vs k-factor
fixpar.C1.NLoS.ds_kf = 0;               % delay spread vs k-factor
fixpar.C1.NLoS.sf_kf = 0;               % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.C1.NLoS.xpr_mu    = 4;          % XPR mean [dB]
fixpar.C1.NLoS.xpr_sigma = 0;          % XPR std  [dB], PK 18.8.2008

%%% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.C1.NLoS.DS_mu      = -7.12;     % delay spread, mean [log10(s)]
fixpar.C1.NLoS.DS_sigma   = 0.33;      % delay spread, std [log10(s)]
fixpar.C1.NLoS.AS_D_mu    = 0.90;      % arrival angle spread, mean [log10(deg)]
fixpar.C1.NLoS.AS_D_sigma = 0.36;      % arrival angle spread, std [log10(deg)]
fixpar.C1.NLoS.AS_A_mu    = 1.65;      % departure angle spread, mean [log10(deg)]
fixpar.C1.NLoS.AS_A_sigma = 0.25;      % departure angle spread, std [log10(deg)]
fixpar.C1.NLoS.SF_sigma   = 8;         % shadowing std [dB] (zero mean)
fixpar.C1.NLoS.KF_mu      = 0;         % k-factor, dummy value
fixpar.C1.NLoS.KF_sigma   = 0;         % k-factor, dummy value

%%% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.C1.NLoS.DS_lambda   = 40;       % [m], delay spread
fixpar.C1.NLoS.AS_D_lambda = 30;       % [m], departure azimuth spread
fixpar.C1.NLoS.AS_A_lambda = 30;       % [m], arrival azimuth spread
fixpar.C1.NLoS.SF_lambda   = 50;       % [m], shadowing
fixpar.C1.NLoS.KF_lambda   = NaN;      % [m], k-factor 

% Path loss, Note! see the path loss equation...
fixpar.C1.NLoS.PL_A = NaN;              % path loss exponent, [d<d_bp d>d_bp]
fixpar.C1.NLoS.PL_B = NaN;              % path loss intercept, [d<d_bp d>d_bp]
fixpar.C1.NLoS.PL_C = NaN;              % path loss frequency dependence factor [dB]
fixpar.C1.NLoS.PL_range = [10 5000];    % applicability range [m], (min max)
fixpar.C1.NLoS.BuildingHeight = 10;     % average building height
%%


%% UMa C2, LoS    
% Fixed scenario specific parameters [1, Table 2]
fixpar.C2.LoS.NumClusters = 12;       % Number of ZDSC
fixpar.C2.LoS.r_DS   = 2.5;           % delays spread proportionality factor
fixpar.C2.LoS.PerClusterAS_D = 5;     % Per cluster FS angle spread [deg]
fixpar.C2.LoS.PerClusterAS_A = 11;    % Per cluster MS angle spread [deg]
fixpar.C2.LoS.LNS_ksi = 3;            % ZDSC LNS ksi [dB], per cluster shadowing

% Cross correlation coefficients [1, Table 2]
fixpar.C2.LoS.asD_ds = 0.4;           % departure AS vs delay spread
fixpar.C2.LoS.asA_ds = 0.8;           % arrival AS vs delay spread
fixpar.C2.LoS.asA_sf = -0.5;          % arrival AS vs shadowing std
fixpar.C2.LoS.asD_sf = -0.5;          % departure AS vs shadowing std
fixpar.C2.LoS.ds_sf  = -0.4;          % delay spread vs shadowing std
fixpar.C2.LoS.asD_asA = 0;            % departure AS vs arrival AS
fixpar.C2.LoS.asD_kf = 0;             % departure AS vs k-factor
fixpar.C2.LoS.asA_kf = -0.2;          % arrival AS vs k-factor, PK 18.8.08
fixpar.C2.LoS.ds_kf = -0.4;           % delay spread vs k-factor
fixpar.C2.LoS.sf_kf = 0;              % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.C2.LoS.xpr_mu    = 8;          % XPR mean [dB]
fixpar.C2.LoS.xpr_sigma = 0;          % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.C2.LoS.DS_mu      = -7.03;     % delay spread, mean [s-dB]
fixpar.C2.LoS.DS_sigma   = 0.66;      % delay spread, std [s-dB]
fixpar.C2.LoS.AS_D_mu    = 1.15;      % arrival angle spread, mean [deg-dB]
fixpar.C2.LoS.AS_D_sigma = 0.28;      % arrival angle spread, std [deg-dB]
fixpar.C2.LoS.AS_A_mu    = 1.81;      % departure angle spread, mean [deg-dB]
fixpar.C2.LoS.AS_A_sigma = 0.20;      % departure angle spread, std [deg-dB]
fixpar.C2.LoS.SF_sigma   = 4;         % shadowing std [dB] (zero mean)
fixpar.C2.LoS.KF_mu      = 9;         % k-factor, dummy value
fixpar.C2.LoS.KF_sigma   = 3.5;       % k-factor, dummy value

% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.C2.LoS.DS_lambda   = 30;        % [m], delay spread
fixpar.C2.LoS.AS_D_lambda = 18;        % [m], departure azimuth spread
fixpar.C2.LoS.AS_A_lambda = 15;        % [m], arrival azimuth spread
fixpar.C2.LoS.SF_lambda   = 37;        % [m], shadowing
fixpar.C2.LoS.KF_lambda   = 12;        % [m], k-factor 

% Path loss, Note! see the path loss equation...
fixpar.C2.LoS.PL_A = [22 40];        % path loss exponent [dB]
fixpar.C2.LoS.PL_B = [28 7.8];       % path loss intercept [dB]
fixpar.C2.LoS.PL_C = [20 2];         % path loss frequency dependence factor
fixpar.C2.LoS.PL_range = [10 5000];  % applicability range [m]

%%


%% UMa C2, NLoS    
% Fixed scenario specific parameters [1, Table 2]
fixpar.C2.NLoS.NumClusters = 20;       % Number of ZDSC
fixpar.C2.NLoS.r_DS   = 2.3;           % delays spread proportionality factor
fixpar.C2.NLoS.PerClusterAS_D = 2;     % Per cluster FS angle spread [deg]
fixpar.C2.NLoS.PerClusterAS_A = 15;    % Per cluster MS angle spread [deg]
fixpar.C2.NLoS.LNS_ksi = 3;            % ZDSC LNS ksi [dB], per cluster shadowing

% Cross correlation coefficients [1, Table 2]
fixpar.C2.NLoS.asD_ds = 0.4;           % departure AS vs delay spread
fixpar.C2.NLoS.asA_ds = 0.6;           % arrival AS vs delay spread
fixpar.C2.NLoS.asA_sf = 0;             % arrival AS vs shadowing std
fixpar.C2.NLoS.asD_sf = -0.6;          % departure AS vs shadowing std
fixpar.C2.NLoS.ds_sf  = -0.4;          % delay spread vs shadowing std
fixpar.C2.NLoS.asD_asA = 0.4;          % departure AS vs arrival AS
fixpar.C2.NLoS.asD_kf = 0;             % departure AS vs k-factor
fixpar.C2.NLoS.asA_kf = 0;             % arrival AS vs k-factor
fixpar.C2.NLoS.ds_kf = 0;              % delay spread vs k-factor
fixpar.C2.NLoS.sf_kf = 0;              % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.C2.NLoS.xpr_mu    = 7;          % XPR mean [dB]
fixpar.C2.NLoS.xpr_sigma = 0;          % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.C2.NLoS.DS_mu      = -6.44;     % delay spread, mean [s-dB]
fixpar.C2.NLoS.DS_sigma   = 0.39;      % delay spread, std [s-dB]
fixpar.C2.NLoS.AS_D_mu    = 1.41;      % arrival angle spread, mean [deg-dB]
fixpar.C2.NLoS.AS_D_sigma = 0.28;      % arrival angle spread, std [deg-dB]
fixpar.C2.NLoS.AS_A_mu    = 1.87;      % departure angle spread, mean [deg-dB]
fixpar.C2.NLoS.AS_A_sigma = 0.11;      % departure angle spread, std [deg-dB]
fixpar.C2.NLoS.SF_sigma   = 6;         % shadowing std [dB] (zero mean)
fixpar.C2.NLoS.KF_mu      = 0;         % k-factor, dummy value
fixpar.C2.NLoS.KF_sigma   = 0;         % k-factor, dummy value

% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.C2.NLoS.DS_lambda   = 40;       % [m], delay spread
fixpar.C2.NLoS.AS_D_lambda = 50;       % [m], departure azimuth spread
fixpar.C2.NLoS.AS_A_lambda = 50;       % [m], arrival azimuth spread
fixpar.C2.NLoS.SF_lambda   = 50;       % [m], shadowing
fixpar.C2.NLoS.KF_lambda   = NaN;      % [m], k-factor 

% Path loss, Note! see the path loss equation...
fixpar.C2.NLoS.PL_A = NaN;            % path loss exponent [dB]
fixpar.C2.NLoS.PL_B = NaN;            % path loss intercept [dB]
fixpar.C2.NLoS.PL_C = NaN;            % path loss frequency dependence factor
fixpar.C2.NLoS.PL_range = [50 5000];  % applicability range [m]
fixpar.C2.NLoS.BuildingHeight = 20;   % average building height

%%


%% UMi B4, NLoS = O-to-I
% Fixed scenario specific parameters
fixpar.B4.NLoS.NumClusters = 12;        % Number of ZDSC    [1, Table 2]
fixpar.B4.NLoS.r_DS   = 2.2;            % delays spread proportionality factor
fixpar.B4.NLoS.PerClusterAS_D = 5;      % Per cluster FS angle spread [deg] [1, Table 2]
fixpar.B4.NLoS.PerClusterAS_A = 8;      % Per cluster MS angle spread [deg] [1, Table 2]
fixpar.B4.NLoS.LNS_ksi = 4;             % ZDSC LNS ksi [dB], per cluster shadowing [1, Table 2]

% Cross correlation coefficients [1, Table 2]
fixpar.B4.NLoS.asD_ds = 0.4;            % departure AS vs delay spread
fixpar.B4.NLoS.asA_ds = 0.4;            % arrival AS vs delay spread
fixpar.B4.NLoS.asA_sf = 0.0;            % arrival AS vs shadowing std
fixpar.B4.NLoS.asD_sf = 0.2;            % departure AS vs shadowing std
fixpar.B4.NLoS.ds_sf  = -0.5;           % delay spread vs shadowing std
fixpar.B4.NLoS.asD_asA = 0;             % departure AS vs arrival AS
fixpar.B4.NLoS.asD_kf = 0;              % departure AS vs k-factor
fixpar.B4.NLoS.asA_kf = 0;              % arrival AS vs k-factor
fixpar.B4.NLoS.ds_kf = 0;               % delay spread vs k-factor
fixpar.B4.NLoS.sf_kf = 0;               % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.B4.NLoS.xpr_mu    = 9;           % XPR mean [dB]
fixpar.B4.NLoS.xpr_sigma = 0;           % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.B4.NLoS.DS_mu      = -6.62;      % delay spread, mean [log10(s)]
fixpar.B4.NLoS.DS_sigma   = 0.32;       % delay spread, std [log10(s)]
fixpar.B4.NLoS.AS_D_mu    = 1.25;       % arrival angle spread, mean [log10(deg)]
fixpar.B4.NLoS.AS_D_sigma = 0.42;       % arrival angle spread, std [log10(deg)]
fixpar.B4.NLoS.AS_A_mu    = 1.76;       % departure angle spread, mean [log10(deg)]
fixpar.B4.NLoS.AS_A_sigma = 0.16;       % departure angle spread, std [log10(deg)]
fixpar.B4.NLoS.SF_sigma   = 7;          % shadowing std [dB] (zero mean)
fixpar.B4.NLoS.KF_mu      = 0;          % k-factor, dummy value
fixpar.B4.NLoS.KF_sigma   = 0;          % k-factor, dummy value

% "Decorrelation distances" [1, Table 2]
fixpar.B4.NLoS.DS_lambda   = 10;         % [m], delay spread
fixpar.B4.NLoS.AS_D_lambda = 11;         % [m], departure azimuth spread
fixpar.B4.NLoS.AS_A_lambda = 17;         % [m], arrival azimuth spread
fixpar.B4.NLoS.SF_lambda   = 7;          % [m], shadowing
fixpar.B4.NLoS.KF_lambda   = NaN;        % [m], k-factor 

% Path loss, Note! see the path loss equation...
fixpar.B4.NLoS.PL_A = NaN;     % path loss exponent 
fixpar.B4.NLoS.PL_B = NaN;     % path loss intercept
fixpar.B4.NLoS.PL_C = 23;               % path loss frequency dependence factor
fixpar.B4.NLoS.PL_range = [3 1000];     % applicability range [m], (min max)
%%


%% RMa D1, LoS
%%% Fixed scenario specific parameters [1, Table 2]
fixpar.D1.LoS.NumClusters = 11;       % Number of ZDSC
fixpar.D1.LoS.r_DS   = 3.8;           % delays spread proportionality factor
fixpar.D1.LoS.PerClusterAS_D = 2;     % Per cluster FS angle spread [deg]
fixpar.D1.LoS.PerClusterAS_A = 3;     % Per cluster MS angle spread [deg]
fixpar.D1.LoS.LNS_ksi = 3;            % ZDSC LNS ksi [dB], per cluster shadowing

% Cross correlation coefficients [1, Table 2]
fixpar.D1.LoS.asD_ds = 0;           % departure AS vs delay spread
fixpar.D1.LoS.asA_ds = 0;           % arrival AS vs delay spread
fixpar.D1.LoS.asA_sf = 0;           % arrival AS vs shadowing std
fixpar.D1.LoS.asD_sf = 0;           % departure AS vs shadowing std
fixpar.D1.LoS.ds_sf  = -0.5;        % delay spread vs shadowing std
fixpar.D1.LoS.asD_asA = 0;          % departure AS vs arrival AS
fixpar.D1.LoS.asD_kf = 0;           % departure AS vs k-factor
fixpar.D1.LoS.asA_kf = 0;           % arrival AS vs k-factor
fixpar.D1.LoS.ds_kf = 0;            % delay spread vs k-factor
fixpar.D1.LoS.sf_kf = 0;            % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.D1.LoS.xpr_mu    = 12;         % XPR mean [dB]
fixpar.D1.LoS.xpr_sigma = 0;          % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.D1.LoS.DS_mu      = -7.49;     % delay spread, mean [s-dB]
fixpar.D1.LoS.DS_sigma   = 0.55;      % delay spread, std [s-dB]
fixpar.D1.LoS.AS_D_mu    = 0.90;      % arrival angle spread, mean [deg-dB]
fixpar.D1.LoS.AS_D_sigma = 0.38;      % arrival angle spread, std [deg-dB]
fixpar.D1.LoS.AS_A_mu    = 1.52;      % departure angle spread, mean [deg-dB]
fixpar.D1.LoS.AS_A_sigma = 0.24;      % departure angle spread, std [deg-dB]
fixpar.D1.LoS.SF_sigma   = 4;         % shadowing std [dB] (zero mean)
fixpar.D1.LoS.KF_mu = 7;              % k-factor mean [dB]
fixpar.D1.LoS.KF_sigma = 4;           % k-factor std [dB]

% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.D1.LoS.DS_lambda   = 50;         % [m], delay spread
fixpar.D1.LoS.AS_D_lambda = 25;         % [m], departure azimuth spread
fixpar.D1.LoS.AS_A_lambda = 35;         % [m], arrival azimuth spread
fixpar.D1.LoS.SF_lambda   = 37;         % [m], shadowing
fixpar.D1.LoS.KF_lambda   = 40;         % [m], k-factor 

% Path loss, Note! see the path loss equation...
fixpar.D1.LoS.PL_A = NaN;               % path loss exponent [dB]
fixpar.D1.LoS.PL_B = NaN;               % path loss intercept [dB]
fixpar.D1.LoS.PL_C = NaN;               % path loss frequency dependence factor
fixpar.D1.LoS.PL_range = [10 10000];    % applicability range [m]
fixpar.D1.LoS.BuildingHeight = 5;       % average building height



%% RMa D1, NLoS
% Fixed scenario specific parameters [1, Table 2]
fixpar.D1.NLoS.NumClusters = 10;        % Number of ZDSC
fixpar.D1.NLoS.r_DS   = 1.7;            % delays spread proportionality factor
fixpar.D1.NLoS.PerClusterAS_D = 2;      % Per cluster FS angle spread [deg]
fixpar.D1.NLoS.PerClusterAS_A = 3;      % Per cluster MS angle spread [deg]
fixpar.D1.NLoS.LNS_ksi = 3;             % ZDSC LNS ksi [dB], per cluster shadowing

% Cross correlation coefficients [1, Table 2]
fixpar.D1.NLoS.asD_ds = -0.4;           % departure AS vs delay spread
fixpar.D1.NLoS.asA_ds = 0;              % arrival AS vs delay spread
fixpar.D1.NLoS.asA_sf = 0;              % arrival AS vs shadowing std
fixpar.D1.NLoS.asD_sf = 0.6;            % departure AS vs shadowing std
fixpar.D1.NLoS.ds_sf  = -0.5;           % delay spread vs shadowing std
fixpar.D1.NLoS.asD_asA = 0;             % departure AS vs arrival AS
fixpar.D1.NLoS.asD_kf = 0;              % departure AS vs k-factor
fixpar.D1.NLoS.asA_kf = 0;              % arrival AS vs k-factor
fixpar.D1.NLoS.ds_kf = 0;               % delay spread vs k-factor
fixpar.D1.NLoS.sf_kf = 0;               % shadowing std vs k-factor

% Polarisation parameters [1, Table 2]
fixpar.D1.NLoS.xpr_mu    = 7;           % XPR mean [dB]
fixpar.D1.NLoS.xpr_sigma = 0;           % XPR std  [dB], PK 18.8.2008

% Dispersion parameters [1, Table 2]
% Log-normal distributions
fixpar.D1.NLoS.DS_mu      = -7.43;      % delay spread, mean [s-dB]
fixpar.D1.NLoS.DS_sigma   = 0.48;       % delay spread, std [s-dB]
fixpar.D1.NLoS.AS_D_mu    = 0.95;       % arrival angle spread, mean [deg-dB]
fixpar.D1.NLoS.AS_D_sigma = 0.45;       % arrival angle spread, std [deg-dB]
fixpar.D1.NLoS.AS_A_mu    = 1.52;       % departure angle spread, mean [deg-dB]
fixpar.D1.NLoS.AS_A_sigma = 0.13;       % departure angle spread, std [deg-dB]
fixpar.D1.NLoS.SF_sigma   = 8;          % shadowing std [dB] (zero mean)
fixpar.D1.NLoS.KF_mu      = 0;          % k-factor, dummy value
fixpar.D1.NLoS.KF_sigma   = 0;          % k-factor, dummy value

%%% Decorrelation distances: lambda parameters [1, Table 2]
fixpar.D1.NLoS.DS_lambda   = 36;        % [m], delay spread
fixpar.D1.NLoS.AS_D_lambda = 30;        % [m], departure azimuth spread
fixpar.D1.NLoS.AS_A_lambda = 40;        % [m], arrival azimuth spread
fixpar.D1.NLoS.SF_lambda   = 120;       % [m], shadowing
fixpar.D1.NLoS.KF_lambda   = NaN;       % [m], k-factor 

% Path loss, Note! see the path loss equation...
fixpar.D1.NLoS.PL_A = NaN;              % path loss exponent [dB]
fixpar.D1.NLoS.PL_B = NaN;              % path loss intercept [dB]
fixpar.D1.NLoS.PL_C = NaN;              % path loss frequency dependence factor
fixpar.D1.NLoS.PL_range = [10 5000];    % applicability range [m]
fixpar.D1.NLoS.BuildingHeight = 5;      % average building height

%%

