function [offset_matrix_AoD, offset_matrix_AoA] = offset_matrix_generation(offset,iterpar)    

%OFFSET_MATRIX_GENERATION Reshapes the offsetmatrix in correct form
%
% Author: Mikko Alatossava (CWC/UOULU)
%
% Date: 14.2.2007
%

N(1) = iterpar.N(1);
N(2) = iterpar.N(2);
N_max = max(N(1),N(2));
NumLinks = length(iterpar.UserIndeces);
LoSConnectionLinks = iterpar.LoSConnectionLinks;
NLoSConnectionLinks = iterpar.NLoSConnectionLinks;
NumLoSConnectionLinks = iterpar.NumLoSConnectionLinks;
NumNLoSConnectionLinks = iterpar.NumNLoSConnectionLinks;

%reshape to cover all the clusters in columns and all the M subpaths in rows
offset_matrix = [offset;-offset];
offset_matrix = repmat(offset_matrix(:),1,N_max*NumLinks); %[M x N*NumLinks] matrix

%find the columns of LoS and NLoS paths
if LoSConnectionLinks
    LoS_Columns = repmat([(LoSConnectionLinks-1)*N_max].',1,N_max)+repmat([1:N_max],NumLoSConnectionLinks,1);
else
   LoS_Columns = [];
end
if NLoSConnectionLinks
    NLoS_Columns = repmat([(NLoSConnectionLinks-1)*N_max].',1,N_max)+repmat([1:N_max],NumNLoSConnectionLinks,1);
else
    NLoS_Columns = [];
end

%find columns in the matrix that are not used, NaN columns
if (N(1) < N(2)) & ~isempty(LoS_Columns)
    NaN_Columns = LoS_Columns(:,end-(N(2)-N(1))+1:end);
elseif (N(2) < N(1)) & ~isempty(NLoS_Columns)
    NaN_Columns = NLoS_Columns(:,end-(N(1)-N(2))+1:end);
else
    NaN_Columns = [];
end

LoS_Columns = reshape(LoS_Columns.',1,[]);
NLoS_Columns = reshape(NLoS_Columns.',1,[]);

%different PerClusterAS_D/A multiplication for LoS/NLoA and AoD/AoA
if LoSConnectionLinks
    offset_matrix_AoD(:,LoS_Columns) = iterpar.LoS.PerClusterAS_D*offset_matrix(:,LoS_Columns);  %[1, eq. 4.14]
    offset_matrix_AoA(:,LoS_Columns) = iterpar.LoS.PerClusterAS_A*offset_matrix(:,LoS_Columns);  %[1, eq. 4.14]
end
if NLoSConnectionLinks
    offset_matrix_AoD(:,NLoS_Columns) = iterpar.NLoS.PerClusterAS_D*offset_matrix(:,NLoS_Columns);  %[1, eq. 4.14]
    offset_matrix_AoA(:,NLoS_Columns) = iterpar.NLoS.PerClusterAS_A*offset_matrix(:,NLoS_Columns);  %[1, eq. 4.14]
end

%set NaN_Columsn to NaN, not used in the calculation
offset_matrix_AoD(:,NaN_Columns) = NaN;
offset_matrix_AoA(:,NaN_Columns) = NaN;
