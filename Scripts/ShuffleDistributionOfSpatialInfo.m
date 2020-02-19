function SpatialInfoMat=ShuffleDistributionOfSpatialInfo(SpatialOccupancy,NumberOfEventsVec,NumOfShuffles)
% SpatialOccupancy=SpatialOccupancy,NumberOfEventsVec;
% NumberOfEventsVec=[50 70];
% NumOfShuffles=1000;

if size(SpatialOccupancy,1)>1
    SpatialOccupancy=SpatialOccupancy';
end

if size(SpatialOccupancy,1)>1 & size(SpatialOccupancy,2)>1
    disp('Problematic size of SpatialOccupancy')
end
        
OccupancyP=SpatialOccupancy;
OccupancyCumP=cumsum(SpatialOccupancy)/sum(SpatialOccupancy);

SpatialInfoMat=nan(NumOfShuffles,length(NumberOfEventsVec));
for runNumOfEvents=1:length(NumberOfEventsVec)
    NumOfEvents=NumberOfEventsVec(runNumOfEvents);
    if NumOfEvents==1
        SpatialInfoMat(:,runNumOfEvents)=NaN;       
    else
    ShufflesMat=histc(rand(NumOfEvents,NumOfShuffles),[0 OccupancyCumP])';
    
    ShufflesMat=ShufflesMat(:,1:end-1);
    MapsOfShuffles=ShufflesMat./repmat(SpatialOccupancy,[NumOfShuffles 1]);
    
    MapsOfShuffles(isnan(MapsOfShuffles))=0;
    MapsOfShuffles=MapsOfShuffles+eps;
    
    mean_r=NumOfEvents/sum(SpatialOccupancy);
    SpatialInfoVec=sum(log2(MapsOfShuffles/mean_r).*(MapsOfShuffles/mean_r).*repmat(OccupancyP,[NumOfShuffles 1]),2);
    SpatialInfoMat(:,runNumOfEvents)=SpatialInfoVec;
    end
end
