clear, clc
%NofMs = 10; NofBs = 6; SectPerBs = 3; K = 10;
d = 100;
%linkpar = layout2link(layoutparset(NofMs,NofBs,SectPerBs,K,d)); % K = # of active links
linkpar = linkparset(24);
antpar = antparset;
wimpar = wimparset;
wimpar.CenterFrequency = 5.25e9;
wimpar.NumTimeSamples = 1e4;
wimpar.NumBsElements = 2;
wimpar.NumMsElements = 2;
wimpar.PolarisedArrays = 'no';
wimpar.PathLossModelUsed = 'yes';
wimpar.ShadowingModelUsed = 'yes';
%linkpar.ScenarioVector = 3*ones(1,10);
linkpar.ScenarioVector = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6];
%linkpar.PropagConditionVector = zeros(1,24);
linkpar.PropagConditionVector = [1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0 0 0 1 0 0 1];
[h,delay,full] = wim(wimpar,linkpar,antpar);
full
