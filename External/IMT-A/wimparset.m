function wimpar=wimparset(varargin)
%WIMPARSET Model parameter configuration for WIMi
%   WIMPAR=WIMPARSET sets default parameters for the input struct WIMPAR 
%   (see WIM). 
%
%   WIMPARSET parameters [ {default} ]:
%
%   NumBsElements           - Number of BS array antenna elements [ {2} ]
%   NumMsElements           - Number of MS array antenna elements [ {2} ]
%   range                   - if Scenario='B5b', the path-loss ranges 1, 2 and 3 are defined, see [1]
%   end_time                - if Scenario='B5x', Observation end time for B5 - time points are taken as [0,T]
%   SampleDensity           - number of time samples per half wavelength [ {2} ]
%   NumTimeSamples          - number of time samples [ {100} ]
%   UniformTimeSampling     - Use same time sampling grid for all links [ yes | {no} ] 
%   IntraClusterDsUsed      - Two strongest clusters are divided into three subclusters [ {yes} | no ] 
%   NumSubPathsPerPath      - number of subpaths per path [ {20} ] (cannot be changed)
%   FixedPdpUsed            - nonrandom path delays and powers [ yes | {no}]
%   FixedAnglesUsed         - nonrandom AoD/AoAs [ yes | {no} ]
%   PolarisedArrays         - usage of dual polarised arrays [ yes | {no} ]
%   CenterFrequency         - carrier frequency in Herz [ {5.25e9} ]
%   DelaySamplingInterval   - delay sampling grid [ {5e-9} ]
%   PathLossModelUsed       - usage of path loss model [ yes | {no} ]
%   ShadowingModelUsed      - usage of shadow fading model [ yes | {no} ]
%   PathLossModel           - path loss model function name [ {pathloss} ]
%   AnsiC_core              - use optimized computation [ yes | {no} ]
%   LookUpTable             - look up EXP(j*THETA) from a table [{0}]
%   RandomSeed              - sets random seed [ {[empty]} ]
%   UseManualPropCondition  - whether to use manual propagation condition (los/nlos) setting or not. 
%                             If not, the propagation condition is drawn from probabilities.  
%
%   Notes about parameters:
%   - The number of BS and MS elements is normally extracted from ANTPAR.
%     The values of NumBsElements and NumMsElements are used only if a single
%     scalar is given as the antenna field pattern in ANTPAR (see ANTPARSET).
%   - For successful Doppler analysis, one should select SampleDensity > 1.
%     The time sample interval is calculated from CenterFrequency and
%     MsVelocity (see LINKPARSET) according to wavelength/(MsVelocity*SampleDensity).
%     The calculated time sample interval for each link is included in the optional 
%     output argument of WIM. 
%   - If UniformTimeSampling is 'yes' all links will be sampled at
%     simultaneous time instants. In this case, the time sample interval is
%     the same for all links it is calculated by replacing MsVelocity with
%     MAX(MsVelocity), where the maximum is over all links. 
%   - Number of rays is fixed to 20. This is because the AoD/AoAs for
%     subpaths in WIM have fixed angle spread. see [1, Table 4-1].
%   - If FixedPdpUsed='yes', the delays and powers of paths are taken from
%     a table [1, Table 6-1..26].  
%   - If FixedAnglesUsed='yes', the AoD/AoAs are taken from a table 
%     [1, Table 6-1..26]. Random pairing of AoDs and AoAs is not used.
%   - If PolarisedArrays='yes', single channel coefficient of impulse
%     response turns to 2x2 coefficient matrix, with elements [VV VH;HV HH].
%     Where V stands for vertical polarisation and H for horizontal.
%   - CenterFrequency affects path loss and time sampling interval.
%   - DelaySamplingInterval determines the sampling grid in delay domain.
%     All path delays are rounded to the nearest grid point. It can also 
%     be set to zero. 
%   - When PathLossModelUsed is 'no' the path losses are still computed for
%     each link but they are not multiplied into the channel matrices. If
%     ShadowingModelUsed is also 'no', each channel matrix element has unit
%     mean power (summed over delay domain). In other words,
%     MEAN(MEAN(ABS(SUM(H,3)).^2,4),5) is a matrix of (approximately) ones 
%     when isotropic unit-gain antennas are used. Exception: with
%     'polarized' option (and default antennas) the mean power is two.
%   - Path loss model is implemented in a separate function, whose name is
%     defined in PathLossModel. For syntax, see PATHLOSS. 
%   - The C-function must be compiled before usage. For more information 
%     of the ANSI-C core function, see SCM_MEX_CORE. 
%   - The LookUpTable parameter defines the number of points used in the
%     cosine look-up table; a power-of-2 should be given. The look-up table 
%     is used only in the ANSI-C optimized core function. Value 0 indicates
%     that look-up table is not used. Value -1 uses the default number of 
%     points, which is 2^14=16384. Since a large part of computation in WIM 
%     involves repeated evaluation of a complex exponential, the look-up 
%     table can speed up computation on certain platforms and C compilers.
%   - Even fixing the random seed may not result in fully repeatable
%     simulations due to differences in e.g. MATLAB versions. 
%
%   Ref. [1]: D1.1.2 V1.0, "WINNER II channel models"
%        [2]: 3GPP TR 25.996 v6.1.0 (2003-09)
%
%   See also WIM, LINKPARSET, LAYOUTPARSET, ANTPARSET.

%   Authors: Jari Salo (HUT), Pekka Kyösti (EBIT), Daniela Laselva (EBIT), 
%   Giovanni Del Galdo (TUI), Marko Milojevic (TUI), Christian Schneider (TUI)
%   Lassi Hentilä (EBIT), Mikko Alatossva (CWC/UOULU)


if length(varargin)>0
    error('No such functionality yet. Try ''wimpar=wimparset'' instead.')
end

% Set the default values
wimpar=struct(  'NumBsElements',2,...                   
                'NumMsElements',2,...                   
                'range',1,...  
                'end_time',1,...                        % Observation end time for B5 - time points are taken as:  wimpar.TimeVector=linspace(0,wimpar.end_time,T);
                'SampleDensity', 2,...                  % in samples/half-wavelength
                'NumTimeSamples',100,...         
                'UniformTimeSampling','no',... 
                'IntraClusterDsUsed','yes',...          % Two strongest clusters are divided into three subclusters
                'NumSubPathsPerPath',20,...             % only value supported is 20.
                'FixedPdpUsed','no',...                 % Use fixed delays and path powers
                'FixedAnglesUsed','no',...              % Use fixed AoD/AoAs
                'PolarisedArrays','no',...              % use polarised arrays
                'TimeEvolution','no',...                % use of time evolution option
                'CenterFrequency',5.25e9,...            % in Herz
                'DelaySamplingInterval',5e-9,...        
                'PathLossModelUsed','no',...            
                'ShadowingModelUsed','no',...           
                'PathLossModel','pathloss',...
                'AnsiC_core','no',...                   
                'LookUpTable',0,...                     % number of points in Ansi-C core look-up table for cosine, 0 if not used
                'RandomSeed',[],...                     % if empty, seed is not set. 
                'UseManualPropCondition','yes');
