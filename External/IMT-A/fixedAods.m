%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed AoDs for different scenarios
% Needed when wimpar.FixedAnglesUsed='yes'
function [aods_los,aod_clusterAS_los, ...
          aods_nlos,aod_clusterAS_nlos]=fixedAods(wimpar,iterpar)
    
Scenario = iterpar.Scenario;

switch Scenario

     case {'A2'}
        
        aods_nlos = [-34 34 60 24 -1 -31 55 69 34 -73 -55 56 -47 -84 87 71 -99 -99 86 ];
        aod_clusterAS_nlos = 5;      % Cluster ASD [deg], [1, table 7-3]
        
        aods_los = [0 79 -84 -97 -102 80 69 104 98 98 -113 -79 -109 89 95];
        aod_clusterAS_los = 5;      
        
    
    case {'B1'}
        
        aods_los = [0 31 35 -38 45 -46 46 -43 -47 51 -62 -58];
        aod_clusterAS_los = 3;      % Cluster ASD [deg], [1, table 7-4]

        aods_nlos = [31 10 26 -27 -3 29 -33 37 -29 27 -35 -48 -50 -44 46 -45 -59 60 56];
        aod_clusterAS_nlos = 10;      % Cluster ASD [deg], [1, table 7-5]

       
     case {'B4'}%need parameters

        aods_nlos = [0 32 -21 37 -43 28 -49 -34 -49 43 49 51];
        aod_clusterAS_nlos = 5;      % Cluster ASD [deg], [1, table 6-3]
        
        aods_los = NaN;
        aod_clusterAS_los = NaN;      
        
        
    case {'C1'}

        aods_los = [0 -29 -32 -31 31 29 -33 35 -30 35 -32 35 33 34 35];
        aod_clusterAS_los = 5;      % Cluster ASD [deg], [1, table 7-10]

        aods_nlos = [0 13 -15 -8 12 -17 12 -8 -10 -13 12 8 14 22];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 7-11]

    case {'C2'}
        aods_los = [0 36 -28 -25 42 -40 -28 -35 -47 -39 47 -43];
        aod_clusterAS_los = 5; 
        
        aods_nlos = [6 44 2 -34 26 -41 -17 -33 24 -34 -38 44 53 54 53 52 57 53 -54 -60];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 7-12]
        

    case {'D1'}

        aods_los = [0 24 23 24 -25 -23 -25 27 27 25 -26];
        aod_clusterAS_los = 2;      % Cluster ASD [deg], [1, table 7-14]

        aods_nlos = [-12 16 0 9 -9 13 -9 13 15 18];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 7-15]
      

end % switch