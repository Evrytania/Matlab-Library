%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed AoAs for different scenarios
% Needed when wimpar.FixedAnglesUsed='yes'
function [aoas_los,aoa_clusterAS_los, ...
          aoas_nlos,aoa_clusterAS_nlos] = fixedAoas(wimpar,iterpar)

Scenario = iterpar.Scenario;

switch Scenario

    case {'A2'}
        aoas_nlos = [50 50 -83 -30 -3 44 -63 -105 -61 110 -85 -85 58 121 -109 108 139 -137 124];
        aoa_clusterAS_nlos = 11;      % Cluster ASA [deg], [1, table 7-3]
        
        aoas_los = [0 85 73 -82 -72 -64 -73 -83 96 -91 -93 74 -87 -77 -97];
        aoa_clusterAS_los = 8; 
        
    case {'B1'}
        aoas_los = [0 -108 153 -111 179 -150 -163 166 14 -4 43 44];
        aoa_clusterAS_los = 17;      % Cluster ASA [deg], [1, table 7-4]

        aoas_nlos = [-83 54 37 96 3 -86 84 -112 -78 65 -88 -130 153 106 147 118 -178 -163 163];
        aoa_clusterAS_nlos = 22;      % Cluster ASA [deg], [1, table 7-5]

    case {'B4'}%need parameters

        aoas_nlos = [0 102 -66 -119 139 91 157 -111 157 138 158 165];
        aoa_clusterAS_nlos = 8;      % Cluster ASA [deg], [1, table 7-9]

        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 

    case {'C1'}
        aoas_los = [0 -144 -159 155 156 -146 168 -176 149 -176 -159 -176 -165 -171 177];
        aoa_clusterAS_los = 5;      % Cluster ASA [deg], [1, table 7-10]

        aoas_nlos = [0 -71 -84 46 -66 -97 -66 -46 -56 73 70 -46 -80 123];
        aoa_clusterAS_nlos = 10;      % Cluster ASA [deg], [1, table 7-11]

    case {'C2'}
        aoas_los = [0 143 -156 133 180 19 161 3 40 -37 -36 -48];
        aoa_clusterAS_los = 11; 
        
        aoas_nlos = [29 -98 8 -114 70 107 59 -103 73 -111 -112 122 129 153 -145 -157 -178 -114 -160 -175];
        aoa_clusterAS_nlos = 15;      % Cluster ASA [deg], [1, table 7-12]
     
    case {'D1'}        

        aoas_los = [0 99 95 99 -106 96 -103 113 110 106 -108];
        aoa_clusterAS_los = 3;      % Cluster ASA [deg], [1, table 7-14]

        aoas_nlos = [-46 58 0 34 34 49 -34 49 56 67];
        aoa_clusterAS_nlos = 3;      % Cluster ASA [deg], [1, table 7-15]


end % switch