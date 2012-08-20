function output = ConstantCMapping(input)
%CONSTANTCMAPPING Determines the constant C used in aod/aoa calculations 
%   from the number of cluster.

%   Authors: Mikko Alatossava (CWC/UOULU)
%   Created: 12.2.2007

switch (input)
    case {4}
        output = 0.779;
    case {5}
        output = 0.860;
    case {8}
        output = 1.018;
    case {10}
        output = 1.090;
    case {11}
        output = 1.123;
    case {12}
        output = 1.146;
    case {14}
        output = 1.190;
    case {15}
        output = 1.211;
    case {16}
        output = 1.226;
    case {20}
        output = 1.289;
    case {22}
        output = 1.300; %tdb, here so that B3 (fixed pdp on) works without problems
    otherwise
        error('Number of Clusters must be among [4 5 8 10 11 12 14 15 16 20]!')                    
end