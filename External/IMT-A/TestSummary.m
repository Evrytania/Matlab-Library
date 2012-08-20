% #1 To test RMa Scenario,LOS propagation condtion when MsBsDistance is longer than breakpoint
% distance.
wimpar=wimparset;linkpar=linkparset;linkpar.MsBsDistance=8000;linkpar.ScenarioVector=6;
linkpar.PropagConditionVector=1;antpar=antparset;wim(wimpar,linkpar,antpar)
