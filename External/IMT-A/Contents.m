% IMT.EVAL channel model
%   based on Winner Phase II channel model
% Version 0.2, September 22, 2008
%
% Channel model functions
%   wim               - WINNER Phase II channel model (D1.1.2) 
%   wimparset         - Model parameter configuration for WIM
%   linkparset        - Link parameter configuration for WIM
%   layoutparset      - Layout parameter configuration for WIM (optional)
%   antparset         - Antenna parameter configuration for WIM
%   pathloss          - Pathloss models for 2GHz and 5GHz 
%
% WINNER -specific functions
%   scenpartables     - Set WIM parameters for WINNER scenarios
%   
% Miscellaneous functions
%   cas               - Circular angle spread (3GPP TR 25.996)
%   ds                - RMS delay spread 
%   dipole            - Field pattern of half wavelength dipole
%   NTlayout          - Visualisation of network layout
%
% Utility functions
%   interp_gain       - Antenna field pattern interpolation
%   interp_gain_c     - Antenna field pattern interpolation (requires GSL)
%   wim_core          - Channel coefficient computation for a geometric channel model
%   scm_mex_core      - WIM_CORE written in ANSI-C 
%   generate_bulk_par - Generation of WIM bulk parameters
%   layout2link       - Computes and converts layout to link parameters
%   ScenarioMapping   - Maps scenario names (A1 etc.) to number indices
%   struct_generation - Assistant function
%   offset_matrix_generation - Assistant function