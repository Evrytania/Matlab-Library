function output = ScenarioMapping(input)
%SCENARIOMAPPING Scenario is given as numerical vector with one element representing 
%   each link. This function maps the element to the corresponding letters as 
%   1=InH (A2), 2=UMi (B1), 3=SMa (C1), 4=UMa (C2), 5=UMi O-2-I (B4) and 6=RMa (D1)

%   Authors: Mikko Alatossava (CWC/UOULU)
%   Created: 12.2.2007
%   Modifications: Updated to IMT.EVAL scenarios      19.9.2008 PekKy

switch (input)
    case {1}
        output = 'A2'; %InH
    case {2}
        output = 'B1'; %UMi
    case {3}
        output = 'C1'; %SMa
    case {4}
        output = 'C2'; %UMa
    case {5}
        output = 'B4'; %UMiO2I
    case {6}
        output = 'D1'; %RMa
end