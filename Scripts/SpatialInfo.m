function SpatialInfoVec=SpatialInfo(SpatialOccupancy,FiringRateMaps)
% This function computes the spatial information for each cell.

% Inputs:
% 1. SpatialOccupancy - a vector of size M (number of spatial bins) with fraction of time spent in each spatial bin (should some up to 1) 
% 2. FiringRateMaps - a matrix of size NxM, where N is the number of cells and M is the number of spatial bins
% each element represents the activity rate of a given cell in a given spatial bin

% Outputs:
% 1. SpatialInfoVec - a vector of size N with the spatial information computed for each cell

OccupancyP=SpatialOccupancy;
  
FiringRateMaps(isnan(FiringRateMaps))=0;
FiringRateMaps=FiringRateMaps+10^(-10);

if size(FiringRateMaps,2)==1
    mean_r=FiringRateMaps'*OccupancyP;
    NormalizedFiringRateMaps=FiringRateMaps'./repmat(mean_r,[1 size(FiringRateMaps,1)]);
    SpatialInfoVec=sum(log2(NormalizedFiringRateMaps).* NormalizedFiringRateMaps.*repmat(OccupancyP',[size(FiringRateMaps,2) 1]),2);
else
    mean_r=FiringRateMaps*OccupancyP;
    NormalizedFiringRateMaps=FiringRateMaps./repmat(mean_r,[1 size(FiringRateMaps,2)]);
    SpatialInfoVec=sum(log2(NormalizedFiringRateMaps).* NormalizedFiringRateMaps.*repmat(OccupancyP',[size(FiringRateMaps,1) 1]),2);
end

