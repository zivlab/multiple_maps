function SpatialInfoVec=SpatialInfo(bins_occupancy,FiringRateMaps)


OccupancyP=bins_occupancy;
  
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

