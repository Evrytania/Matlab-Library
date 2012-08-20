function bulk_parameters = struct_generation(bulk_parameters, ...
                                            bulk_parameters_iter, ...
                                            wimpar,linkpar,iterpar, ...
                                            state) 

% STRUCT_GENERATION Generates, updates and refines struct bulk_parameters
%
% Author: Mikko Alatossava (CWC/UOULU)
%
% Date: 14.2.2007
%
% Revision history after 12.2.2007:
%
% XPRv and XPRh unified     10.10.2007 HentLas


%%
MsBsDistance          = linkpar.MsBsDistance;

if strcmp(upper(state),'INITIALIZATION') 
    clear bulk_parameters
    %reserve space in bulkp_params-struct for all scenario dependent-struct
    %max number of clusters is 20. At the end the N-columns that have NaNs are discarded
    MaxClusterValue = 20;
    bulk_parameters.delays = NaN*ones(length(MsBsDistance),MaxClusterValue);
    bulk_parameters.path_powers = NaN*ones(length(MsBsDistance),MaxClusterValue);
    bulk_parameters.aods = NaN*ones(length(MsBsDistance),MaxClusterValue,wimpar.NumSubPathsPerPath);
    bulk_parameters.aoas = NaN*ones(length(MsBsDistance),MaxClusterValue,wimpar.NumSubPathsPerPath);
    bulk_parameters.path_losses = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.MsBsDistance = NaN*ones(1,length(MsBsDistance));
    bulk_parameters.shadow_fading = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.o2v_shadow_fading = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.sigmas = NaN*ones(length(MsBsDistance),5);
    bulk_parameters.propag_condition = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.LoSO2ILinks = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.LoSO2VLinks = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.NLoSO2ILinks = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.NLoSO2VLinks = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.Kcluster = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.Phi_LOS = NaN*ones(length(MsBsDistance),1);
    bulk_parameters.scatterer_freq = NaN*ones(length(MsBsDistance),MaxClusterValue,wimpar.NumSubPathsPerPath);
    if strcmp(wimpar.PolarisedArrays,'no')
        bulk_parameters.subpath_phases = NaN*ones(length(MsBsDistance),MaxClusterValue,wimpar.NumSubPathsPerPath);
    elseif strcmp(wimpar.PolarisedArrays,'yes')
          bulk_parameters.subpath_phases = NaN*ones(length(MsBsDistance),4,MaxClusterValue,wimpar.NumSubPathsPerPath);
          bulk_parameters.xpr = NaN*ones(length(MsBsDistance),MaxClusterValue,wimpar.NumSubPathsPerPath);
    end

elseif strcmp(upper(state),'ITERATION')
    bulk_parameters.delays(bulk_parameters_iter.user_indeces, 1:size(bulk_parameters_iter.delays,2))...
        = bulk_parameters_iter.delays;
    bulk_parameters.path_powers(bulk_parameters_iter.user_indeces, 1:size(bulk_parameters_iter.path_powers,2))... 
        = bulk_parameters_iter.path_powers;
    bulk_parameters.aods(bulk_parameters_iter.user_indeces,1:size(bulk_parameters_iter.aods,2),:) = bulk_parameters_iter.aods;
    bulk_parameters.aoas(bulk_parameters_iter.user_indeces,1:size(bulk_parameters_iter.aoas,2),:) = bulk_parameters_iter.aoas;
    bulk_parameters.path_losses(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.path_losses;
    bulk_parameters.MsBsDistance(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.MsBsDistance;
    bulk_parameters.shadow_fading(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.shadow_fading;
    bulk_parameters.o2v_shadow_fading(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.o2v_shadow_fading;
    bulk_parameters.propag_condition(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.propag_condition;
    bulk_parameters.LoS02ILinks(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.LoS02ILinks(bulk_parameters_iter.user_indeces);
    bulk_parameters.LoS02VLinks(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.LoS02VLinks(bulk_parameters_iter.user_indeces);
    bulk_parameters.NLoS02ILinks(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.NLoS02ILinks(bulk_parameters_iter.user_indeces);
    bulk_parameters.NLoS02VLinks(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.NLoS02VLinks(bulk_parameters_iter.user_indeces);
    bulk_parameters.Kcluster(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.Kcluster;
    bulk_parameters.Phi_LOS(bulk_parameters_iter.user_indeces) = bulk_parameters_iter.Phi_LOS;
    if strcmp(iterpar.Scenario(1:2),'B5')
        bulk_parameters.scatterer_freq(bulk_parameters_iter.user_indeces,1:size(bulk_parameters_iter.scatterer_freq,2),:) = bulk_parameters_iter.scatterer_freq; 
    else
        bulk_parameters.sigmas(bulk_parameters_iter.user_indeces,:) = bulk_parameters_iter.sigmas;
    end
    if strcmp(wimpar.PolarisedArrays,'no')
        bulk_parameters.subpath_phases(bulk_parameters_iter.user_indeces,1:size(bulk_parameters_iter.subpath_phases,2),:) = bulk_parameters_iter.subpath_phases;
    elseif strcmp(wimpar.PolarisedArrays,'yes')
        bulk_parameters.subpath_phases(bulk_parameters_iter.user_indeces,:,1:size(bulk_parameters_iter.subpath_phases,3),:) = bulk_parameters_iter.subpath_phases;
        bulk_parameters.xpr(bulk_parameters_iter.user_indeces,1:size(bulk_parameters_iter.xpr,2),:) = bulk_parameters_iter.xpr;
    end
    
elseif strcmp(upper(state),'REFINEMENT')
    while all(isnan(bulk_parameters.delays(:,end))) & all(isnan(bulk_parameters.aods(:,end,:)))
        bulk_parameters.delays = bulk_parameters.delays(:,1:end-1);
        bulk_parameters.path_powers = bulk_parameters.path_powers(:,1:end-1);
        bulk_parameters.aods = bulk_parameters.aods(:,1:end-1,:);
        bulk_parameters.aoas = bulk_parameters.aoas(:,1:end-1,:);
        if isfield(bulk_parameters,'scatterer_freq')
            bulk_parameters.scatterer_freq = bulk_parameters.scatterer_freq(:,1:end-1,:); 
        end
        if strcmp(wimpar.PolarisedArrays,'no')
            bulk_parameters.subpath_phases = bulk_parameters.subpath_phases(:,1:end-1,:);
        elseif strcmp(wimpar.PolarisedArrays,'yes')
            bulk_parameters.subpath_phases = bulk_parameters.subpath_phases(:,:,1:end-1,:);
            bulk_parameters.xpr = bulk_parameters.xpr(:,1:end-1,:);
        end
    end
end