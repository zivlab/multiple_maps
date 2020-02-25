function SpatialInfoMat=ShuffleDistributionOfSpatialInfo(SpatialOccupancy,NumberOfEventsVec,NumOfShuffles)
% This function computes the spatial information for each cell.

% Inputs:
% 1. SpatialOccupancy - a vector of size M (number of spatial bins) with fraction of time spent in each spatial bin (should some up to 1) 
% 2. NumberOfEventsVec - a vector of size K with all the unique numbers of events for
% each cell - each unique event number requires a separate shuffle. 
% 3. NumOfShuffles - the number of shuffles run for each cell (cells with
% the same number of events will have the same shuffles).

% Outputs:
% 1. SpatialInfoMat - a matrix of size NumOfShuffles X K. All shuffled
% values of spatial information are provided for each unique number of events.
% The decision of which cells are significantly tuned to position according
% to a defined p-value takes place outside this function.

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
