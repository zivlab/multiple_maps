%% Setting the data paths and parameters:
clc
clear all
close all

%%
% Choose mouse and environment:
chosen_mouse='C7_M4'; 
chosen_environment='environment B'; %environment B;
data_pathway='D:\dev\multiple_maps\multiple_maps\Data\';
mkdir([data_pathway chosen_mouse '\' chosen_environment],'multiple_maps_results')
mkdir([data_pathway chosen_mouse '\' chosen_environment '\multiple_maps_results'],'figures')
figures_directory=[data_pathway chosen_mouse '\' chosen_environment '\multiple_maps_results\figures'];

% loading the map numbers per trial:
maps_per_trial=xlsread([data_pathway chosen_mouse '\' chosen_environment '\maps per trial.xlsx']);
map_number_per_trial=maps_per_trial(2:end,2:end);

velocity_thresh=1; % in cm/sec units - previously used: 3;
track_length=96;
num_bins=24;
bins_to_use=3:22;
relevant_trials=2:6;
fs=20;
velocity_kernel_length=0.4*fs;
bin_size=track_length/num_bins;
num_sessions=size(map_number_per_trial,1);
num_trials=length(relevant_trials);

%% Loading the data
% parameters:

% Loading the data:
behavior=cell(1,num_sessions);
frame_log=cell(1,num_sessions);
neuronal_activity=cell(1,num_sessions);
for n=1:num_sessions
    neuronal_activity{n}=load([data_pathway chosen_mouse '\' chosen_environment '\finalEventsMat_' num2str(n) '.mat']);    
    frame_log{n}=xlsread([data_pathway chosen_mouse '\' chosen_environment '\frameLog_' num2str(n) '.csv']);
    behavior{n}=load([data_pathway  chosen_mouse '\' chosen_environment '\behavior_' num2str(n)]);
end
if num_sessions>1
    registration_results=load(fullfile([data_pathway chosen_mouse '\' chosen_environment],'cell_registration.mat'));
    cell_to_index_map=registration_results.cell_to_index_map;
    number_of_neurons=size(cell_to_index_map,1);
else
    number_of_neurons=size(neuronal_activity{1}.allEventsMat,2);
end


%%
% binning and smoothing the behavioral data:
velocity_smoothing_kernel=1/velocity_kernel_length*ones(1,velocity_kernel_length);
for n=1:num_sessions
    for k=1:num_trials
        if isfield(behavior{n}.my_mvmt{relevant_trials(k)},'velocity')
            temp_velocity=behavior{n}.my_mvmt{relevant_trials(k)}.velocity;
        else
            temp_velocity=behavior{n}.my_mvmt{relevant_trials(k)}.velocity_VER2;
        end
        smoothed_velocity=conv(temp_velocity,velocity_smoothing_kernel,'same')./conv(ones(length(temp_velocity),1),velocity_smoothing_kernel,'same');
        if isfield(behavior{n}.my_mvmt{relevant_trials(k)},'velocity')
            behavior{n}.my_mvmt{relevant_trials(k)}.velocity=smoothed_velocity;
        else
            behavior{n}.my_mvmt{relevant_trials(k)}.velocity_VER2=smoothed_velocity;
        end
    end
end

% reverse bins for L-shape:
reverse_bins=1:13;
for n=1:num_sessions
    for k=1:num_trials
        if ~isfield(behavior{n}.my_mvmt{relevant_trials(k)},'velocity')
            temp_locations=behavior{n}.my_mvmt{relevant_trials(k)}.binned;
            for m=1:length(reverse_bins)
                behavior{n}.my_mvmt{relevant_trials(k)}.binned(temp_locations==m)=max(reverse_bins)+1-m;
            end
        end
    end
end

% Creating a concataned location vector for two directions of movement:
all_locations=cell(num_sessions,num_trials);
all_velocities=cell(num_sessions,num_trials);
all_events=cell(num_sessions,num_trials);
neurons_per_session=zeros(1,num_sessions);
for n=1:num_sessions
    for k=1:num_trials
        minimal_position=0;
        maximal_position=1;
        if isfield(behavior{n}.my_mvmt{relevant_trials(k)},'position')
            minimal_position=min([minimal_position, min(behavior{n}.my_mvmt{relevant_trials(k)}.position(2:end))]);
            maximal_position=max([maximal_position, max(behavior{n}.my_mvmt{relevant_trials(k)}.position(2:end))]);
        else
            minimal_position=min([minimal_position, min(behavior{n}.my_mvmt{relevant_trials(k)}.rough_linearized_location(2:end))]);
            maximal_position=max([maximal_position, max(behavior{n}.my_mvmt{relevant_trials(k)}.rough_linearized_location(2:end))]);
        end
        position_range=[minimal_position maximal_position];
        if isfield(behavior{n}.my_mvmt{relevant_trials(k)},'position')
            all_locations{n,k}=ceil((behavior{n}.my_mvmt{relevant_trials(k)}.position((2:end))-position_range(1))./(position_range(2)-position_range(1))*num_bins);
            all_locations{n,k}(all_locations{n,k}==0)=1;
        else
            all_locations{n,k}=ceil((behavior{n}.my_mvmt{relevant_trials(k)}.rough_linearized_location((2:end))-position_range(1))./(position_range(2)-position_range(1))*num_bins);
            all_locations{n,k}(all_locations{n,k}==0)=1;
        end
        if isfield(behavior{n}.my_mvmt{relevant_trials(k)},'velocity')
            all_velocities{n,k}=behavior{n}.my_mvmt{relevant_trials(k)}.velocity(2:end);
        else
            all_velocities{n,k}=behavior{n}.my_mvmt{relevant_trials(k)}.velocity_VER2(2:end);
        end
        this_event_index=frame_log{n}(relevant_trials(k),3):frame_log{n}(relevant_trials(k),4);
        if num_sessions==1
            all_events{n,k}=neuronal_activity{n}.allEventsMat(this_event_index,:);
        else % for registered data:
            if isfield(neuronal_activity{n},'finalSpikesMat')
                this_session_events=full(neuronal_activity{n}.finalSpikesMat(:,this_event_index))';
                this_session_events(this_session_events<6.5)=0;
            else
                this_session_events=neuronal_activity{n}.allEventsMat(this_event_index,:);
            end
            registered_events=zeros(size(this_session_events,1),size(cell_to_index_map,1));
            num_cells=size(this_session_events,2);
            neurons_per_session(n)=num_cells;
            for m=1:num_cells
                this_registered_index=find(cell_to_index_map(:,n)==m);
                if ~isempty(this_registered_index)
                    registered_events(:,this_registered_index)=this_session_events(:,m);
                end
            end
            all_events{n,k}=registered_events;
        end
    end
end

% matching the numbers of frames per trial for events, locations and
% velocities:
number_of_frames_per_trial=zeros(num_sessions,num_trials);
binary_events=cell(num_sessions,num_trials);
number_of_calcium_frames=zeros(num_sessions,num_trials);
number_of_behavioral_frames=zeros(num_sessions,num_trials);
for n=1:num_sessions
    for k=1:num_trials
        number_of_frames_per_trial(n,k)=min([size(all_velocities{n,k},1),size(all_locations{n,k},1),size(all_events{n,k},1)]);
        number_of_calcium_frames(n,k)=size(all_events{n,k},1);
        number_of_behavioral_frames(n,k)=size(all_locations{n,k},1);
        all_events{n,k}=all_events{n,k}(1:number_of_frames_per_trial(n,k),:);
        binary_events{n,k}=all_events{n,k}>0;
        all_velocities{n,k}=all_velocities{n,k}(1:number_of_frames_per_trial(n,k));
        all_locations{n,k}=all_locations{n,k}(1:number_of_frames_per_trial(n,k));
    end
end

% Constructing the prior distributions:
prior_per_trial=zeros(num_sessions,num_trials,2,length(bins_to_use));
prior_per_session=zeros(num_sessions,2,length(bins_to_use));
number_of_frames_per_trial_running=zeros(num_sessions,2,num_trials);
for n=1:num_sessions
    for k=1:num_trials
        positive_veloctity=find(all_velocities{n,k}>=velocity_thresh);
        negative_veloctity=find(all_velocities{n,k}<-velocity_thresh);
        for m=1:length(bins_to_use)
            prior_per_trial(n,k,1,m)=prior_per_trial(n,k,1,m)+sum(all_locations{n,k}(positive_veloctity)==bins_to_use(m));
            prior_per_trial(n,k,2,m)=prior_per_trial(n,k,2,m)+sum(all_locations{n,k}(negative_veloctity)==bins_to_use(m));
        end
        number_of_frames_per_trial_running(n,1,k)=length(positive_veloctity);
        number_of_frames_per_trial_running(n,2,k)=length(negative_veloctity);
        prior_per_trial(n,k,1,:)=prior_per_trial(n,k,1,:)./sum(prior_per_trial(n,k,1,:));
        prior_per_trial(n,k,2,:)=prior_per_trial(n,k,2,:)./sum(prior_per_trial(n,k,2,:));
    end
    prior_per_session(n,:,:)=squeeze(sum(squeeze(prior_per_trial(n,:,:,:)).*repmat(squeeze(number_of_frames_per_trial_running(n,:,:))',1,1,length(bins_to_use)),1))./repmat(sum(squeeze(number_of_frames_per_trial_running(n,:,:)),2),1,length(bins_to_use));
end

%% Computing rate correlations and overlap:
overlap_events_threshold=5;
number_of_maps=max(max(map_number_per_trial));

rate_correlations=zeros(num_sessions*num_trials,num_sessions*num_trials);
shuffle_rate_correlations=zeros(num_sessions*num_trials,num_sessions*num_trials);
overlap=zeros(num_sessions*num_trials,num_sessions*num_trials);
shuffle_overlap=zeros(num_sessions*num_trials,num_sessions*num_trials);
all_rate_vectors=zeros(num_sessions*num_trials,number_of_neurons);
number_events_trials=zeros(num_sessions,num_trials,number_of_neurons);
for n=1:num_sessions*num_trials
    first_session_index=ceil(n/num_trials);
    first_trial_index=mod(n,num_trials);
    if first_trial_index==0
        first_trial_index=num_trials;
    end
    all_rate_vectors(n,:)=sum(all_events{first_session_index,first_trial_index});
    number_events_trials(first_session_index,first_trial_index,:)=sum(all_events{first_session_index,first_trial_index}>0);
    for k=1:num_sessions*num_trials
        second_session_index=ceil(k/num_trials);
        second_trial_index=mod(k,num_trials);
        if second_trial_index==0
            second_trial_index=num_trials;
        end
        first_session_events=all_events{first_session_index,first_trial_index};
        first_session_binary_activity=sum(first_session_events>0);
        second_session_events=all_events{second_session_index,second_trial_index};
        second_session_binary_activity=sum(second_session_events>0);
        [~,random_order]=sort(rand(1,size(first_session_binary_activity,2)));        
        shuffle_first_session_binary_activity=first_session_binary_activity(random_order);
        rate_correlations(n,k)=corr(first_session_binary_activity',second_session_binary_activity');
        shuffle_rate_correlations(n,k)=corr(shuffle_first_session_binary_activity',second_session_binary_activity');
        overlap(n,k)=sum(first_session_binary_activity>=overlap_events_threshold & second_session_binary_activity>=overlap_events_threshold)/sum(first_session_binary_activity>=overlap_events_threshold | second_session_binary_activity>=overlap_events_threshold);
        shuffle_overlap(n,k)=sum(shuffle_first_session_binary_activity>=overlap_events_threshold & second_session_binary_activity>=overlap_events_threshold)/sum(shuffle_first_session_binary_activity>=overlap_events_threshold | second_session_binary_activity>=overlap_events_threshold);
    end
end
concatanated_events=cell(1,num_sessions);
concatanated_binary_events=cell(1,num_sessions);
concatanated_locations=cell(1,num_sessions);
concatanated_velocities=cell(1,num_sessions);
for k=1:num_sessions
    concatanated_events{k}=[];
    concatanated_locations{k}=[];
    concatanated_velocities{k}=[];
    concatanated_binary_events{k}=[];
    for n=1:num_trials
        concatanated_events{k}=[concatanated_events{k}; all_events{k,n}];
        concatanated_locations{k}=[concatanated_locations{k}; all_locations{k,n}];
        concatanated_velocities{k}=[concatanated_velocities{k}; all_velocities{k,n}];
        concatanated_binary_events{k}=[concatanated_binary_events{k}; binary_events{k,n}];
    end
end

% PCA for rate vectors:
reshaped_maps_vector=reshape(map_number_per_trial',num_trials*num_sessions,1);
temp_color_map=colormap('jet');
close;
color_map=temp_color_map(4:size(temp_color_map,1)/num_sessions:end,:);
[~,score,~]=pca(all_rate_vectors);
color_vector=[1 0 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 0 1 1 ; 1 1 0];

% long-term_dynamics:
within_maps_rate_dynamics=cell(1,num_sessions);
across_maps_rate_dynamics=cell(1,num_sessions);
shuffle_maps_rate_dynamics=cell(1,num_sessions);
all_maps_rate_dynamics=cell(1,num_sessions);
within_maps_overlap_dynamics=cell(1,num_sessions);
across_maps_overlap_dynamics=cell(1,num_sessions);
shuffle_maps_overlap_dynamics=cell(1,num_sessions);
all_maps_overlap_dynamics=cell(1,num_sessions);
for n=1:num_sessions*num_trials
    first_session=ceil(n/num_trials);
    first_map=reshaped_maps_vector(n);
    for k=n+1:num_sessions*num_trials
        second_session=ceil(k/num_trials);
        second_map=reshaped_maps_vector(k);
        all_maps_rate_dynamics{second_session-first_session+1}=[all_maps_rate_dynamics{second_session-first_session+1},rate_correlations(n,k)];
        all_maps_overlap_dynamics{second_session-first_session+1}=[all_maps_overlap_dynamics{second_session-first_session+1},overlap(n,k)];
        if first_map==second_map
            within_maps_rate_dynamics{second_session-first_session+1}=[within_maps_rate_dynamics{second_session-first_session+1},rate_correlations(n,k)];
            within_maps_overlap_dynamics{second_session-first_session+1}=[within_maps_overlap_dynamics{second_session-first_session+1},overlap(n,k)];
        else
            across_maps_rate_dynamics{second_session-first_session+1}=[across_maps_rate_dynamics{second_session-first_session+1},rate_correlations(n,k)];
            shuffle_maps_rate_dynamics{second_session-first_session+1}=[shuffle_maps_rate_dynamics{second_session-first_session+1},shuffle_rate_correlations(n,k)];
            across_maps_overlap_dynamics{second_session-first_session+1}=[across_maps_overlap_dynamics{second_session-first_session+1},overlap(n,k)];
            shuffle_maps_overlap_dynamics{second_session-first_session+1}=[shuffle_maps_overlap_dynamics{second_session-first_session+1},shuffle_overlap(n,k)];
        end
    end
end
mean_within_maps_rate_dynamics=zeros(1,num_sessions);
mean_across_maps_rate_dynamics=zeros(1,num_sessions);
mean_shuffle_rate_dynamics=zeros(1,num_sessions);
mean_all_maps_rate_dynamics=zeros(1,num_sessions);
std_within_maps_rate_dynamics=zeros(1,num_sessions);
std_across_maps_rate_dynamics=zeros(1,num_sessions);
std_shuffle_rate_dynamics=zeros(1,num_sessions);
std_all_maps_rate_dynamics=zeros(1,num_sessions);
mean_within_maps_overlap_dynamics=zeros(1,num_sessions);
mean_across_maps_overlap_dynamics=zeros(1,num_sessions);
mean_shuffle_overlap_dynamics=zeros(1,num_sessions);
mean_all_maps_overlap_dynamics=zeros(1,num_sessions);
std_within_maps_overlap_dynamics=zeros(1,num_sessions);
std_across_maps_overlap_dynamics=zeros(1,num_sessions);
std_shuffle_overlap_dynamics=zeros(1,num_sessions);
std_all_maps_overlap_dynamics=zeros(1,num_sessions);
for n=1:num_sessions
    mean_within_maps_rate_dynamics(n)=mean(within_maps_rate_dynamics{n},'omitnan');
    std_within_maps_rate_dynamics(n)=std(within_maps_rate_dynamics{n},'omitnan');
    mean_across_maps_rate_dynamics(n)=mean(across_maps_rate_dynamics{n},'omitnan');
    std_across_maps_rate_dynamics(n)=std(across_maps_rate_dynamics{n},'omitnan');
    mean_all_maps_rate_dynamics(n)=mean(all_maps_rate_dynamics{n},'omitnan');
    std_all_maps_rate_dynamics(n)=std(all_maps_rate_dynamics{n},'omitnan');
    mean_shuffle_rate_dynamics(n)=mean(shuffle_maps_rate_dynamics{n},'omitnan');
    std_shuffle_rate_dynamics(n)=std(shuffle_maps_rate_dynamics{n},'omitnan');
    
    mean_within_maps_overlap_dynamics(n)=mean(within_maps_overlap_dynamics{n},'omitnan');
    std_within_maps_overlap_dynamics(n)=std(within_maps_overlap_dynamics{n},'omitnan');
    mean_across_maps_overlap_dynamics(n)=mean(across_maps_overlap_dynamics{n},'omitnan');
    std_across_maps_overlap_dynamics(n)=std(across_maps_overlap_dynamics{n},'omitnan');
    mean_all_maps_overlap_dynamics(n)=mean(all_maps_overlap_dynamics{n},'omitnan');
    std_all_maps_overlap_dynamics(n)=std(all_maps_overlap_dynamics{n},'omitnan');
    mean_shuffle_overlap_dynamics(n)=mean(shuffle_maps_overlap_dynamics{n},'omitnan');
    std_shuffle_overlap_dynamics(n)=std(shuffle_maps_overlap_dynamics{n},'omitnan');
end

% plotting results:
figure
errorbar(2*(0:num_sessions-1),mean_within_maps_rate_dynamics,std_within_maps_rate_dynamics,'--*b')
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_rate_dynamics,std_across_maps_rate_dynamics,'--*r')
hold on
errorbar(2*(0:num_sessions-1),mean_all_maps_rate_dynamics,std_all_maps_rate_dynamics,'--*k')
errorbar(2*(0:num_sessions-1),mean_shuffle_rate_dynamics,std_shuffle_rate_dynamics,'--*g')
xlim([0 2*(num_sessions-1)])
xlabel('Elapsed time (days)')
ylabel('Rate correlation')
set(gca,'fontsize',16)
legend('Within map','Across maps','All maps','shuffle')
legend('boxoff')

figure
errorbar(2*(0:num_sessions-1),mean_within_maps_overlap_dynamics,std_within_maps_overlap_dynamics,'--*b')
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_overlap_dynamics,std_across_maps_overlap_dynamics,'--*r')
hold on
errorbar(2*(0:num_sessions-1),mean_all_maps_overlap_dynamics,std_all_maps_overlap_dynamics,'--*k')
errorbar(2*(0:num_sessions-1),mean_shuffle_overlap_dynamics,std_shuffle_overlap_dynamics,'--*g')
xlim([0 2*(num_sessions-1)])
xlabel('Elapsed time (days)')
ylabel('Overlap')
set(gca,'fontsize',16)
legend('Within map','Across maps','All maps','shuffle')
legend('boxoff')

figure('units','normalized','outerposition',[0.35 0.2 0.3 0.6])
axes('position',[0.1 0.1 0.75 0.75])
imagesc(rate_correlations)
xlabel('Session number')
ylabel('Session number')
min_correlation=min(min(rate_correlations));
caxis([round(min_correlation*100)/100 1])
set(gca,'xtick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'xticklabel',1:num_sessions)
set(gca,'ytick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'yticklabel',1:num_sessions)
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:num_sessions-1
    hold on
    plot(0.5+[0 num_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','k')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 num_sessions*num_trials],'linewidth',2,'color','k')
end
axes('position',[0.87 0.13 0.03 0.69])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)
savefig(fullfile(figures_directory,'Rate correlations.fig'))
saveas(gcf,fullfile(figures_directory,'Rate correlations'),'png')

%% Example place cells:
% Example cells from mouse C7_M4 environment B:
chosen_session=5;
maps_vector=[1 1 2 1 2];
chosen_cells=[82 62 42 120 4];

% Example cells from mouse C6_M4 environment A:
% chosen_session=1;
% maps_vector=[1 1 1 1 1];
% chosen_cells=[333,22,69,56,16];

number_of_frames_per_map(1)=sum(number_of_frames_per_trial(chosen_session,:));
number_of_frames_per_map(2)=sum(number_of_frames_per_trial(chosen_session,maps_vector==1));
number_of_frames_per_map(3)=sum(number_of_frames_per_trial(chosen_session,maps_vector==2));

events_per_map_per_direction=cell(length(chosen_cells),2,3);
tuning_per_map_per_direction=zeros(length(chosen_cells),2,3,num_bins);
smoothing_kernel=gausswin(1+round(num_bins/4));
chosen_cells_temp=[];
for n=1:length(chosen_cells)
    this_cell_index=find(cell_to_index_map(:,chosen_session)==chosen_cells(n));
    for k=1:num_trials        
        chosen_cell_events=find(all_events{chosen_session,k}(:,this_cell_index)>0);
        event_locations=all_locations{chosen_session,k}(chosen_cell_events);
        event_velocities=all_velocities{chosen_session,k}(chosen_cell_events);
        events_per_map_per_direction{n,1,1}=[events_per_map_per_direction{n,1,1} ; event_locations(event_velocities>0)];
        events_per_map_per_direction{n,2,1}=[events_per_map_per_direction{n,2,1} ; event_locations(event_velocities<=0)];
        if maps_vector(k)==1
            events_per_map_per_direction{n,1,2}=[events_per_map_per_direction{n,1,2} ; event_locations(event_velocities>0)];
            events_per_map_per_direction{n,2,2}=[events_per_map_per_direction{n,2,2} ; event_locations(event_velocities<=0)];
        else
            events_per_map_per_direction{n,1,3}=[events_per_map_per_direction{n,1,3} ; event_locations(event_velocities>0)];
            events_per_map_per_direction{n,2,3}=[events_per_map_per_direction{n,2,3} ; event_locations(event_velocities<=0)];
        end
    end
    temp_tuning=hist(events_per_map_per_direction{n,1,1},1:num_bins);
    tuning_per_map_per_direction(n,1,1,:)=conv(temp_tuning,smoothing_kernel,'same')./conv(ones(1,num_bins),smoothing_kernel,'same');
    temp_tuning=hist(events_per_map_per_direction{n,1,2},1:num_bins);
    tuning_per_map_per_direction(n,1,2,:)=conv(temp_tuning,smoothing_kernel,'same')./conv(ones(1,num_bins),smoothing_kernel,'same');
    temp_tuning=hist(events_per_map_per_direction{n,1,3},1:num_bins);
    tuning_per_map_per_direction(n,1,3,:)=conv(temp_tuning,smoothing_kernel,'same')./conv(ones(1,num_bins),smoothing_kernel,'same');
    temp_tuning=hist(events_per_map_per_direction{n,2,1},1:num_bins);
    tuning_per_map_per_direction(n,2,1,:)=conv(temp_tuning,smoothing_kernel,'same')./conv(ones(1,num_bins),smoothing_kernel,'same');
    temp_tuning=hist(events_per_map_per_direction{n,2,2},1:num_bins);
    tuning_per_map_per_direction(n,2,2,:)=conv(temp_tuning,smoothing_kernel,'same')./conv(ones(1,num_bins),smoothing_kernel,'same');
    temp_tuning=hist(events_per_map_per_direction{n,2,3},1:num_bins);
    tuning_per_map_per_direction(n,2,3,:)=conv(temp_tuning,smoothing_kernel,'same')./conv(ones(1,num_bins),smoothing_kernel,'same');
end

% plotting example cells:
t=0:1/fs:(length(concatanated_locations{chosen_session})-1)/fs;
num_plotted_trials=5;
temporal_indexes=1:num_plotted_trials/num_trials*length(concatanated_locations{chosen_session});
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
for n=1:length(chosen_cells)
    axes('position',[(n-1)*0.185+0.07 0.42 0.165 0.53])
    plot(bin_size*concatanated_locations{chosen_session}(temporal_indexes),max(t)-t(temporal_indexes))
    this_cell_index=find(cell_to_index_map(:,chosen_session)==chosen_cells(n));

    chosen_cell_events=find(concatanated_events{chosen_session}(temporal_indexes,this_cell_index)>0);
    chosen_cell_velocities=concatanated_velocities{chosen_session}(chosen_cell_events);
    
    num_events_chosen_cell=length(chosen_cell_events);
    for k=1:num_events_chosen_cell
        hold on
        if chosen_cell_velocities(k)>0
            plot(bin_size*concatanated_locations{chosen_session}(chosen_cell_events(k)),max(t)-t(chosen_cell_events(k)),'.','markersize',20,'color','r')
        else
            plot(bin_size*concatanated_locations{chosen_session}(chosen_cell_events(k)),max(t)-t(chosen_cell_events(k)),'.','markersize',20,'color',[0 0.8 0])
        end
    end
    if n==1
        ylabel('Time (sec)')
        set(gca,'ytick',0:100:max(t))
        set(gca,'yticklabel',fliplr(0:100:max(t)))
    else
        set(gca,'ytick',[])
    end
    set(gca,'xtick',[])
    ylim(max(t)- [max(t(temporal_indexes)) 0])
    for k=1:num_plotted_trials-1
        hold on
        plot([0 100],max(t)-[t(frame_log{chosen_session}(relevant_trials(k),3)) t(frame_log{chosen_session}(relevant_trials(k),3))],'linewidth',2,'color','k')
    end
    title(['Cell ' num2str(n)],'fontweight','normal')
    set(gca,'fontsize',16)
    
    axes('position',[(n-1)*0.185+0.07 0.1 0.165 0.27])
    plot(-bin_size/2+bin_size*[1:num_bins],fs*(squeeze(tuning_per_map_per_direction(n,1,2,:)+tuning_per_map_per_direction(n,2,2,:)))/number_of_frames_per_map(2),'color','b','linewidth',2)
    hold on
    plot(-bin_size/2+bin_size*[1:num_bins],fs*(squeeze(tuning_per_map_per_direction(n,1,3,:)+tuning_per_map_per_direction(n,2,3,:)))/number_of_frames_per_map(3),'color','m','linewidth',2)
    hold on
    plot(-bin_size/2+bin_size*[1:num_bins],fs*(squeeze(tuning_per_map_per_direction(n,1,1,:)+tuning_per_map_per_direction(n,2,1,:)))/number_of_frames_per_map(1),'--','color','k','linewidth',2)
    max_number_of_events=round(max(squeeze(tuning_per_map_per_direction(n,1,1,:)+tuning_per_map_per_direction(n,2,1,:))));
    if n==1
        xlabel('Position (cm)')
        ylabel('Event rate')
        legend('Map 1','Map 2','Combined')
        legend('boxoff')
        set(gca,'ytick',[0 8])
    else
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    set(gca,'fontsize',16)
    box off
end

%% Place tuning analysis (trial level):
event_threshold=3;

positive_veloctity_trials=cell(num_sessions,num_trials);
negative_veloctity_trials=cell(num_sessions,num_trials);
running_indexes_trials=cell(num_sessions,num_trials);
events_running_trials=cell(num_sessions,num_trials);
binary_events_running_trials=cell(num_sessions,num_trials);
position_running_trials=cell(num_sessions,num_trials);
velocity_running_trials=cell(num_sessions,num_trials);
events_right_trials=cell(num_sessions,num_trials);
binary_events_right_trials=cell(num_sessions,num_trials);
position_right_trials=cell(num_sessions,num_trials);
velocity_right_trials=cell(num_sessions,num_trials);
events_left_trials=cell(num_sessions,num_trials);
binary_events_left_trials=cell(num_sessions,num_trials);
position_left_trials=cell(num_sessions,num_trials);
velocity_left_trials=cell(num_sessions,num_trials);
for k=1:num_sessions
    for n=1:num_trials
        positive_veloctity_trials{k,n}=find(all_velocities{k,n}>=velocity_thresh); % added the location filter after analyzing CA3 data
        negative_veloctity_trials{k,n}=find(all_velocities{k,n}<=-velocity_thresh);
        running_indexes_trials{k,n}=sort(union(positive_veloctity_trials{k,n},negative_veloctity_trials{k,n}));
        if max(running_indexes_trials{k,n})>size(all_events{k,n},1)
            non_existing_indexes=find(running_indexes_trials{k,n}>size(all_events{k,n},1));
            events_running_trials{k,n}=all_events{k,n}(running_indexes_trials{k,n}(1:non_existing_indexes-1),:);
            position_running_trials{k,n}=all_locations{k,n}(running_indexes_trials{k,n}(1:non_existing_indexes-1));
            velocity_running_trials{k,n}=all_velocities{k,n}(running_indexes_trials{k,n}(1:non_existing_indexes-1));
        else
            events_running_trials{k,n}=all_events{k,n}(running_indexes_trials{k,n},:);
            position_running_trials{k,n}=all_locations{k,n}(running_indexes_trials{k,n});
            velocity_running_trials{k,n}=all_velocities{k,n}(running_indexes_trials{k,n});
        end
        binary_events_running_trials{k,n}=events_running_trials{k,n}>0;
        if max(positive_veloctity_trials{k,n})>size(all_events{k,n},1)
            non_existing_indexes=find(positive_veloctity_trials{k,n}>size(all_events{k,n},1));
            events_right_trials{k,n}=all_events{k,n}(positive_veloctity_trials{k,n}(1:non_existing_indexes-1),:);
            position_right_trials{k,n}=all_locations{k,n}(positive_veloctity_trials{k,n}(1:non_existing_indexes-1));
            velocity_right_trials{k,n}=all_velocities{k,n}(positive_veloctity_trials{k,n}(1:non_existing_indexes-1));
        else
            events_right_trials{k,n}=all_events{k,n}(positive_veloctity_trials{k,n},:);
            position_right_trials{k,n}=all_locations{k,n}(positive_veloctity_trials{k,n});
            velocity_right_trials{k,n}=all_velocities{k,n}(positive_veloctity_trials{k,n});
        end
        binary_events_right_trials{k,n}=events_right_trials{k,n}>0;
        
        if max(negative_veloctity_trials{k,n})>size(all_events{k,n},1)
            non_existing_indexes=find(negative_veloctity_trials{k,n}>size(all_events{k,n},1));
            events_left_trials{k,n}=all_events{k,n}(negative_veloctity_trials{k,n}(1:non_existing_indexes-1),:);
            position_left_trials{k,n}=all_locations{k,n}(negative_veloctity_trials{k,n}(1:non_existing_indexes-1));
            velocity_left_trials{k,n}=all_velocities{k,n}(negative_veloctity_trials{k,n}(1:non_existing_indexes-1));
        else
            events_left_trials{k,n}=all_events{k,n}(negative_veloctity_trials{k,n},:);
            position_left_trials{k,n}=all_locations{k,n}(negative_veloctity_trials{k,n});
            velocity_left_trials{k,n}=all_velocities{k,n}(negative_veloctity_trials{k,n});
        end
        binary_events_left_trials{k,n}=events_left_trials{k,n}>0;
    end
end

% computing the rate maps for both directions:
rate_maps_trials=zeros(num_sessions,num_trials,number_of_neurons,length(bins_to_use),2);
for l=1:num_sessions
    for m=1:num_trials
        for k=1:length(bins_to_use)
            this_bins_left_index=find(round(position_left_trials{l,m})==bins_to_use(k));
            this_bins_right_index=find(round(position_right_trials{l,m})==bins_to_use(k));
            for n=1:number_of_neurons
                rate_maps_trials(l,m,n,k,1)=mean(binary_events_right_trials{l,m}(this_bins_right_index,n),'omitnan')';
                rate_maps_trials(l,m,n,k,2)=mean(binary_events_left_trials{l,m}(this_bins_left_index,n),'omitnan')';
            end
        end
    end
end

% smoothing the rate maps:
activity_threshold=0;
smoothing_kernel=gausswin(1+round(length(bins_to_use)/4));
smoothed_rate_maps_trials=zeros(size(rate_maps_trials));
preferred_locations_trials=zeros(num_sessions,num_trials,number_of_neurons,2);
for k=1:num_sessions
    for m=1:num_trials
        for n=1:number_of_neurons
            smoothed_rate_maps_trials(k,m,n,:,1)=conv(squeeze(rate_maps_trials(k,m,n,:,1))',smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_trials(k,m,n,:,2)=conv(squeeze(rate_maps_trials(k,m,n,:,2))',smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            if max(smoothed_rate_maps_trials(k,m,n,:,1))>activity_threshold
                [~,preferred_locations_trials(k,m,n,1)]=max(smoothed_rate_maps_trials(k,m,n,:,1));
            end
            if max(smoothed_rate_maps_trials(k,m,n,:,2))>activity_threshold
                [~,preferred_locations_trials(k,m,n,2)]=max(smoothed_rate_maps_trials(k,m,n,:,2));
            end
        end
    end
end
preferred_locations_trials(preferred_locations_trials==0)=nan;

normalized_rate_maps_trials=zeros(size(rate_maps_trials));
for k=1:num_sessions
    for m=1:num_trials
        for n=1:number_of_neurons
            if max(smoothed_rate_maps_trials(k,m,n,:,1))>0
                normalized_rate_maps_trials(k,m,n,:,1)=smoothed_rate_maps_trials(k,m,n,:,1)/max(smoothed_rate_maps_trials(k,m,n,:,1));
            end
            if max(smoothed_rate_maps_trials(k,m,n,:,2))>0
                normalized_rate_maps_trials(k,m,n,:,2)=smoothed_rate_maps_trials(k,m,n,:,2)/max(smoothed_rate_maps_trials(k,m,n,:,2));
            end
        end
    end
end

ordered_locations_right_trials=cell(num_sessions,num_trials);
ordered_locations_left_trials=cell(num_sessions,num_trials);
not_nan_right_trials=cell(num_sessions,num_trials);
not_nan_left_trials=cell(num_sessions,num_trials);
for k=1:num_sessions
    for m=1:num_trials
        not_nan_right_trials{k,m}=find(~isnan(preferred_locations_trials(k,m,:,1)));
        not_nan_left_trials{k,m}=find(~isnan(preferred_locations_trials(k,m,:,2)));
        [~,ordered_locations_right_trials{k,m}]=sort(preferred_locations_trials(k,m,not_nan_right_trials{k,m},1));
        [~,ordered_locations_left_trials{k,m}]=sort(preferred_locations_trials(k,m,not_nan_left_trials{k,m},2));
    end
end

ordered_normalized_maps_left_trials=cell(num_sessions,num_trials);
for k=1:num_sessions
    for m=1:num_trials
        ordered_normalized_maps_left_trials{k,m}=zeros(length(not_nan_left_trials{k,m}),length(bins_to_use));
        for n=1:length(not_nan_left_trials{k,m})
            ordered_normalized_maps_left_trials{k,m}(n,:)=normalized_rate_maps_trials(k,m,not_nan_left_trials{k,m}(ordered_locations_left_trials{k,m}(n)),:,2);
        end
    end
end
ordered_normalized_maps_right_trials=cell(num_sessions,num_trials);
for k=1:num_sessions
    for m=1:num_trials
        ordered_normalized_maps_right_trials{k,m}=zeros(length(not_nan_right_trials{k,m}),length(bins_to_use));
        for n=1:length(not_nan_right_trials{k,m})
            ordered_normalized_maps_right_trials{k,m}(n,:)=normalized_rate_maps_trials(k,m,not_nan_right_trials{k,m}(ordered_locations_right_trials{k,m}(n)),:,1);
        end
    end
end

% spatial information:
over_threshold_right_trials=cell(num_sessions,num_trials);
rate_maps_right_trials=cell(num_sessions,num_trials);
over_threshold_left_trials=cell(num_sessions,num_trials);
rate_maps_left_trials=cell(num_sessions,num_trials);
active_cells_trials=cell(num_sessions,num_trials);
active_cells_sessions=cell(1,num_sessions);
for k=1:num_sessions
    for m=1:num_trials
        over_threshold_right_trials{k,m}=find(sum(binary_events_right_trials{k,m},1)>=event_threshold);
        rate_maps_right_trials{k,m}=squeeze(rate_maps_trials(k,m,over_threshold_right_trials{k,m},:,1));
        over_threshold_left_trials{k,m}=find(sum(binary_events_left_trials{k,m},1)>=event_threshold);
        rate_maps_left_trials{k,m}=squeeze(rate_maps_trials(k,m,over_threshold_left_trials{k,m},:,2));
        active_cells_trials{k,m}=find(sum(binary_events{k,m},1)>=event_threshold);        
        active_cells_sessions{k}=unique(union(active_cells_sessions{k},active_cells_trials{k,m}));
    end
end

% calculating spatial information:
spatial_information_right_trials=cell(num_sessions,num_trials);
spatial_information_left_trials=cell(num_sessions,num_trials);
spatial_information_trials=cell(num_sessions,num_trials);
for k=1:num_sessions
    for m=1:num_trials
        spatial_information_right_trials{k,m}=SpatialInfo(squeeze(prior_per_trial(k,m,1,:)),squeeze(rate_maps_right_trials{k,m}(:,:)));
        spatial_information_left_trials{k,m}=SpatialInfo(squeeze(prior_per_trial(k,m,2,:)),squeeze(rate_maps_left_trials{k,m}(:,:)));
        spatial_information_trials{k,m}=[spatial_information_left_trials{k,m}',spatial_information_right_trials{k,m}'];
    end
end

% checking for significant place cells:
significant_threshold=0.05;
num_shuffles=1000;
place_significance_left_trials=cell(num_sessions,num_trials);
place_significance_right_trials=cell(num_sessions,num_trials);
shuffle_spatial_information_left_trials=cell(num_sessions,num_trials);
shuffle_spatial_information_right_trials=cell(num_sessions,num_trials);
shuffle_place_significance_left_trials=cell(num_sessions,num_trials);
shuffle_place_significance_right_trials=cell(num_sessions,num_trials);
for k=1:num_sessions
    for m=1:num_trials
        num_events_vec=unique([sum(binary_events_left_trials{k,m}(:,over_threshold_left_trials{k,m}),1), sum(binary_events_right_trials{k,m}(:,over_threshold_right_trials{k,m}),1)]);
        spatial_info_shuffle_right=ShuffleDistributionOfSpatialInfo(squeeze(prior_per_trial(k,m,1,:)),num_events_vec,num_shuffles);
        spatial_info_shuffle_left=ShuffleDistributionOfSpatialInfo(squeeze(prior_per_trial(k,m,2,:)),num_events_vec,num_shuffles);
        
        place_significance_left_trials{k,m}=zeros(1,length(spatial_information_left_trials{k,m}));
        shuffle_spatial_information_left_trials{k,m}=zeros(1,length(spatial_information_left_trials{k,m}));
        for n=1:length(spatial_information_left_trials{k,m})
            this_number_events=sum(binary_events_left_trials{k,m}(:,over_threshold_left_trials{k,m}(n)),1);
            this_number_events_index=find(num_events_vec==this_number_events);
            place_significance_left_trials{k,m}(n)=1-sum(spatial_information_left_trials{k,m}(n)>spatial_info_shuffle_left(:,this_number_events_index))/num_shuffles;
            shuffle_spatial_information_left_trials{k,m}(n)=spatial_info_shuffle_left(randi(num_shuffles,1),this_number_events_index);
            shuffle_place_significance_left_trials{k,m}(n)=1-sum(shuffle_spatial_information_left_trials{k,m}(n)>spatial_info_shuffle_left(:,this_number_events_index))/num_shuffles;
        end
        place_significance_right_trials{k,m}=zeros(1,length(spatial_information_right_trials{k,m}));
        shuffle_spatial_information_right_trials{k,m}=zeros(1,length(spatial_information_right_trials{k,m}));
        for n=1:length(spatial_information_right_trials{k,m})
            this_number_events=sum(binary_events_right_trials{k,m}(:,over_threshold_right_trials{k,m}(n)),1);
            this_number_events_index=find(num_events_vec==this_number_events);
            place_significance_right_trials{k,m}(n)=1-sum(spatial_information_right_trials{k,m}(n)>spatial_info_shuffle_right(:,this_number_events_index))/num_shuffles;
            shuffle_spatial_information_right_trials{k,m}(n)=spatial_info_shuffle_right(randi(num_shuffles,1),this_number_events_index);
            shuffle_place_significance_right_trials{k,m}(n)=1-sum(shuffle_spatial_information_right_trials{k,m}(n)>spatial_info_shuffle_right(:,this_number_events_index))/num_shuffles;
        end
    end
end

preferred_locations_trials_both_sides=zeros(num_sessions,num_trials,2*number_of_neurons);
normalized_rate_maps_trials_both_sides=zeros(num_sessions,num_trials,2*number_of_neurons,length(bins_to_use));
normalized_rate_maps_averaged_both_sides=zeros(num_sessions,num_trials,number_of_neurons,length(bins_to_use));
rate_maps_trials_both_sides=zeros(num_sessions,num_trials,2*number_of_neurons,length(bins_to_use));
smoothed_rate_maps_trials_both_sides=zeros(num_sessions,num_trials,2*number_of_neurons,length(bins_to_use));
for m=1:num_sessions
    preferred_locations_trials_both_sides(m,:,:)=[squeeze(preferred_locations_trials(m,:,:,1)) , squeeze(preferred_locations_trials(m,:,:,2))];
    normalized_rate_maps_trials_both_sides(m,:,:,:)=[squeeze(normalized_rate_maps_trials(m,:,:,:,1)) , squeeze(normalized_rate_maps_trials(m,:,:,:,2))];
    normalized_rate_maps_averaged_both_sides(m,:,:,:)=squeeze(mean([normalized_rate_maps_trials(m,:,:,:,1) ; normalized_rate_maps_trials(m,:,:,:,2)],'omitnan'));
    rate_maps_trials_both_sides(m,:,:,:)=[squeeze(rate_maps_trials(m,:,:,:,1)) , squeeze(rate_maps_trials(m,:,:,:,2))];
    smoothed_rate_maps_trials_both_sides(m,:,:,:)=[squeeze(smoothed_rate_maps_trials(m,:,:,:,1)) , squeeze(smoothed_rate_maps_trials(m,:,:,:,2))];
end

reshaped_smoothed_rate_maps_trials_both_sides=zeros(num_sessions,num_trials,number_of_neurons,length(bins_to_use));
for n=1:num_sessions
    for k=1:num_trials
        for m=1:number_of_neurons
            reshaped_smoothed_rate_maps_trials_both_sides(n,k,m,:)=mean(squeeze([smoothed_rate_maps_trials_both_sides(n,k,m,:), smoothed_rate_maps_trials_both_sides(n,k,number_of_neurons+m,:)]));
        end
    end
end

preferred_locations_averaged_both_sides=zeros(num_sessions,num_trials,number_of_neurons);
for k=1:num_sessions
    for m=1:num_trials
        for n=1:number_of_neurons
            if max(reshaped_smoothed_rate_maps_trials_both_sides(k,m,n,:))>activity_threshold
                [~,preferred_locations_averaged_both_sides(k,m,n)]=max(squeeze(reshaped_smoothed_rate_maps_trials_both_sides(k,m,n,:)));
                normalized_rate_maps_averaged_both_sides(k,m,n,:)=normalized_rate_maps_averaged_both_sides(k,m,n,:)./max(normalized_rate_maps_averaged_both_sides(k,m,n,:));
            end
        end
    end
end

significant_left_trials=cell(num_sessions,num_trials);
significant_right_trials=cell(num_sessions,num_trials);
significant_trials=cell(num_sessions,num_trials);
significant_left_sessions=cell(1,num_sessions);
significant_right_sessions=cell(1,num_sessions);
significant_sessions=cell(1,num_sessions);
fraction_significant_trials=zeros(num_sessions,num_trials);
fraction_significant_sessions=zeros(1,num_sessions);
number_of_significant_cells_trials=zeros(num_sessions,num_trials);
shuffle_number_of_significant_cells_trials=zeros(num_sessions,num_trials);
shuffle_fraction_of_significant_cells_trials=zeros(num_sessions,num_trials);
all_significant=[];
for k=1:num_sessions
    significant_left_sessions{k}=[];
    significant_right_sessions{k}=[];
    for m=1:num_trials
        significant_left_trials{k,m}=over_threshold_left_trials{k,m}(find(place_significance_left_trials{k,m}<significant_threshold));
        significant_right_trials{k,m}=over_threshold_right_trials{k,m}(find(place_significance_right_trials{k,m}<significant_threshold));
        significant_left_sessions{k}=union(significant_left_sessions{k},significant_left_trials{k,m});
        significant_right_sessions{k}=union(significant_right_sessions{k},significant_right_trials{k,m});
        significant_trials{k,m}=[significant_right_trials{k,m} , number_of_neurons+significant_left_trials{k,m}];
        temp_significant_trials=unique([significant_right_trials{k,m} , significant_left_trials{k,m}]);
        fraction_significant_trials(k,m)=length(temp_significant_trials)/length(active_cells_trials{k,m});
        
        number_of_significant_cells_trials(k,m)=length(union(significant_left_trials{k,m},significant_right_trials{k,m}));
        shuffle_significant_left_trials{k,m}=over_threshold_left_trials{k,m}(find(shuffle_place_significance_left_trials{k,m}<significant_threshold));
        shuffle_significant_right_trials{k,m}=over_threshold_right_trials{k,m}(find(shuffle_place_significance_right_trials{k,m}<significant_threshold));
        shuffle_number_of_significant_cells_trials(k,m)=length(union(shuffle_significant_left_trials{k,m},shuffle_significant_right_trials{k,m}));
        shuffle_fraction_of_significant_cells_trials(k,m)=length(union(shuffle_significant_left_trials{k,m},shuffle_significant_right_trials{k,m}))/length(active_cells_trials{k,m});
    end
    significant_sessions{k}=[significant_right_sessions{k} ; number_of_neurons+significant_left_sessions{k}];
    temp_significant_sessions=unique([significant_right_sessions{k} ; significant_left_sessions{k}]);
    if num_sessions>1
        fraction_significant_sessions(k)=length(temp_significant_sessions)/sum(cell_to_index_map(:,k)>0);
    else
        fraction_significant_sessions(k)=length(temp_significant_sessions)/number_of_neurons;
    end
    all_significant=union(all_significant,temp_significant_sessions);
end
fraction_significant=length(all_significant)/(number_of_neurons);

average_spatial_information_trials=zeros(num_sessions,num_trials);
shuffle_spatial_information_trials=zeros(num_sessions,num_trials);
for k=1:num_sessions
    for m=1:num_trials
        temp_significant_left=spatial_information_left_trials{k,m}(find(place_significance_left_trials{k,m}<significant_threshold));
        temp_significant_right=spatial_information_right_trials{k,m}(find(place_significance_right_trials{k,m}<significant_threshold));
        temp_shuffle_left=shuffle_spatial_information_left_trials{k,m}(find(place_significance_left_trials{k,m}<significant_threshold));
        temp_shuffle_right=shuffle_spatial_information_right_trials{k,m}(find(place_significance_right_trials{k,m}<significant_threshold));
        average_spatial_information_trials(k,m)=mean([temp_significant_left; temp_significant_right]);
        shuffle_spatial_information_trials(k,m)=mean([temp_shuffle_left, temp_shuffle_right]);
    end
end

% calculating rate map and PV correlations:
rate_correlations_non_place_cells=ones(num_sessions*num_trials,num_sessions*num_trials);
rate_correlations_place_cells=ones(num_sessions*num_trials,num_sessions*num_trials);
rate_map_correlations=ones(num_sessions*num_trials,num_sessions*num_trials);
positional_shifts=cell(num_sessions*num_trials,num_sessions*num_trials);
PV_correlations=ones(num_sessions*num_trials,num_sessions*num_trials);
PV_correlations_separated=ones(num_sessions*num_trials,num_sessions*num_trials);
shuffle_PV_correlations_separated=ones(num_sessions*num_trials,num_sessions*num_trials);
flipped_PV_correlations=ones(num_sessions*num_trials,num_sessions*num_trials);
PV_correlations_bin_resolution=zeros(num_sessions*num_trials*length(bins_to_use),num_sessions*num_trials*length(bins_to_use));
PV_correlations_bin_resolution_separated=zeros(num_sessions*num_trials*2*length(bins_to_use),num_sessions*num_trials*2*length(bins_to_use));
all_rate_map_correlations=[];
all_PV_correlations=[];
for n=1:num_sessions*num_trials
    for k=n:num_sessions*num_trials
        first_session_index=ceil(n/num_trials);
        first_trial_index=mod(n,num_trials);
        if first_trial_index==0
            first_trial_index=num_trials;
        end
        second_session_index=ceil(k/num_trials);
        second_trial_index=mod(k,num_trials);
        if second_trial_index==0
            second_trial_index=num_trials;
        end
        first_trial_significant=significant_trials{first_session_index,first_trial_index};
        second_trial_significant=significant_trials{second_session_index,second_trial_index};
        this_comparison_significant_cells=union(first_trial_significant,second_trial_significant);
        this_comparison_significant_left=union(significant_left_trials{first_session_index,first_trial_index},significant_left_trials{second_session_index,second_trial_index});
        this_comparison_significant_right=union(significant_right_trials{first_session_index,first_trial_index},significant_right_trials{second_session_index,second_trial_index});
        this_comparison_significant_separated=union(this_comparison_significant_right,this_comparison_significant_left);
        
        this_comparison_union_place_cells=this_comparison_significant_cells;
        this_comparison_union_place_cells(this_comparison_union_place_cells>number_of_neurons)=this_comparison_union_place_cells(this_comparison_union_place_cells>number_of_neurons)-number_of_neurons;
        this_comparison_union_place_cells=unique(this_comparison_union_place_cells);
        
        first_trial_active_cells=active_cells_trials{first_session_index,first_trial_index};
        second_trial_active_cells=active_cells_trials{second_session_index,second_trial_index};
        this_comparison_active_cells=union(first_trial_active_cells,second_trial_active_cells);
        this_comparison_non_place_cells=setdiff(1:number_of_neurons,this_comparison_union_place_cells);
        this_comparison_place_cells=setdiff(1:number_of_neurons,setdiff(this_comparison_active_cells,this_comparison_union_place_cells));
          
        first_session_events=all_events{first_session_index,first_trial_index}(:,this_comparison_non_place_cells);
        first_session_binary_activity=sum(first_session_events>0);
        second_session_events=all_events{second_session_index,second_trial_index}(:,this_comparison_non_place_cells);
        second_session_binary_activity=sum(second_session_events>0);
        rate_correlations_non_place_cells(n,k)=corr(first_session_binary_activity',second_session_binary_activity');
        rate_correlations_non_place_cells(k,n)=rate_correlations_non_place_cells(n,k);
        
        first_session_events=all_events{first_session_index,first_trial_index}(:,this_comparison_place_cells);
        first_session_binary_activity=sum(first_session_events>0);
        second_session_events=all_events{second_session_index,second_trial_index}(:,this_comparison_place_cells);
        second_session_binary_activity=sum(second_session_events>0);
        rate_correlations_place_cells(n,k)=corr(first_session_binary_activity',second_session_binary_activity');
        rate_correlations_place_cells(k,n)=rate_correlations_place_cells(n,k);
        
        all_maps_trial_one=squeeze(smoothed_rate_maps_trials_both_sides(first_session_index,first_trial_index,this_comparison_significant_cells,:));
        all_maps_trial_two=squeeze(smoothed_rate_maps_trials_both_sides(second_session_index,second_trial_index,this_comparison_significant_cells,:));
        all_maps_trial_one_right=squeeze(smoothed_rate_maps_trials(first_session_index,first_trial_index,this_comparison_significant_separated,:,1));
        all_maps_trial_one_left=squeeze(smoothed_rate_maps_trials(first_session_index,first_trial_index,this_comparison_significant_separated,:,2));
        all_maps_trial_two_right=squeeze(smoothed_rate_maps_trials(second_session_index,second_trial_index,this_comparison_significant_separated,:,1));
        all_maps_trial_two_left=squeeze(smoothed_rate_maps_trials(second_session_index,second_trial_index,this_comparison_significant_separated,:,2));
        all_maps_trial_one_separated=[all_maps_trial_one_right , all_maps_trial_one_left];
        all_maps_trial_two_separated=[all_maps_trial_two_right , all_maps_trial_two_left];
        flipped_maps_trial_one_separated=fliplr(all_maps_trial_one_separated);
        
        [~,random_order]=sort(rand(1,size(all_maps_trial_one_separated,1)));
        shuffled_maps_trial_one_separated=all_maps_trial_one_separated(random_order,:);
        
        all_flipped_PV_correlations_temp=corr(flipped_maps_trial_one_separated,all_maps_trial_two_separated);
        all_PV_correlations_separated_temp=corr(all_maps_trial_one_separated,all_maps_trial_two_separated);
        all_shuffle_PV_correlations_separated_temp=corr(shuffled_maps_trial_one_separated,all_maps_trial_two_separated);
        flipped_PV_correlations_temp=zeros(1,2*length(bins_to_use));
        PV_correlations_separated_temp=zeros(1,2*length(bins_to_use));
        shuffle_PV_correlations_separated_temp=zeros(1,2*length(bins_to_use));
        for m=1:2*length(bins_to_use)
            PV_correlations_separated_temp(m)=all_PV_correlations_separated_temp(m,m);
            flipped_PV_correlations_temp(m)=all_flipped_PV_correlations_temp(m,m);
            shuffle_PV_correlations_separated_temp(m)=all_shuffle_PV_correlations_separated_temp(m,m);
        end
        
        all_rate_map_correlations_temp=corr(all_maps_trial_one',all_maps_trial_two');
        rate_map_correlations_temp=zeros(1,length(this_comparison_significant_cells));
        for m=1:length(this_comparison_significant_cells)
            rate_map_correlations_temp(m)=all_rate_map_correlations_temp(m,m);
        end
        all_PV_correlations_temp=corr(all_maps_trial_one,all_maps_trial_two);
        all_PV_correlations_separated=corr(all_maps_trial_one_separated,all_maps_trial_two_separated);
        PV_correlations_temp=zeros(1,length(bins_to_use));
        for m=1:length(bins_to_use)
            PV_correlations_temp(m)=all_PV_correlations_separated(m,m);
        end
        if k~=n
            all_rate_map_correlations=[all_rate_map_correlations , rate_map_correlations_temp(~isnan(rate_map_correlations_temp))];
            all_PV_correlations=[all_PV_correlations , PV_correlations_temp(~isnan(PV_correlations_temp))];
            rate_map_correlations(n,k)=mean(rate_map_correlations_temp,'omitnan');
            rate_map_correlations(k,n)=rate_map_correlations(n,k);
            [~,preferred_postions_map_one]=max(all_maps_trial_one');
            [~,preferred_postions_map_two]=max(all_maps_trial_two');
            positional_shifts{k,n}=preferred_postions_map_one-preferred_postions_map_two;
            positional_shifts{k,n}(max(all_maps_trial_one')<0.02 | max(all_maps_trial_two')<0.02)=[];
            positional_shifts{k,n}(isnan(positional_shifts{k,n}))=[];
            positional_shifts{n,k}=positional_shifts{k,n};
            PV_correlations(n,k)=mean(PV_correlations_temp,'omitnan');
            PV_correlations(k,n)=PV_correlations(n,k);
            PV_correlations_separated(n,k)=mean(PV_correlations_separated_temp,'omitnan');
            PV_correlations_separated(k,n)=PV_correlations_separated(n,k);
            shuffle_PV_correlations_separated(n,k)=mean(shuffle_PV_correlations_separated_temp,'omitnan');
            shuffle_PV_correlations_separated(k,n)=shuffle_PV_correlations_separated(n,k);
            flipped_PV_correlations(n,k)=mean(flipped_PV_correlations_temp,'omitnan');
            flipped_PV_correlations(k,n)=flipped_PV_correlations(n,k);
        end
        PV_correlations_bin_resolution((n-1)*length(bins_to_use)+1:n*length(bins_to_use),(k-1)*length(bins_to_use)+1:k*length(bins_to_use))=all_PV_correlations_temp;
        PV_correlations_bin_resolution((k-1)*length(bins_to_use)+1:k*length(bins_to_use),(n-1)*length(bins_to_use)+1:n*length(bins_to_use))=all_PV_correlations_temp';
        PV_correlations_bin_resolution_separated((n-1)*2*length(bins_to_use)+1:n*2*length(bins_to_use),(k-1)*2*length(bins_to_use)+1:k*2*length(bins_to_use))=all_PV_correlations_separated;
        PV_correlations_bin_resolution_separated((k-1)*2*length(bins_to_use)+1:k*2*length(bins_to_use),(n-1)*2*length(bins_to_use)+1:n*2*length(bins_to_use))=all_PV_correlations_separated';
    end
end

% Eigenvectors of the PV correlation matrix:
[V,D]=eig(PV_correlations);
reduced_PV_correlations=PV_correlations*V;

% plotting the eigenvector analysis of PV correlations:
temp_color_map=colormap('jet');
close;
color_map=temp_color_map(4:size(temp_color_map,1)/num_sessions:end,:);
days_vec=reshape(repmat(1:num_sessions,[num_trials 1]),[1 num_sessions*num_trials]);

figure
for n=1:num_sessions*num_trials
    this_session=ceil(n/num_trials);
    if reshaped_maps_vector(n)>0
        plot(-reduced_PV_correlations(n,num_sessions*num_trials-1),-reduced_PV_correlations(n,num_sessions*num_trials),'.','markersize',20,'color',color_vector(reshaped_maps_vector(n),:))
        hold on
    end
end
xlabel('Component 1')
ylabel('Component 2')
set(gca,'fontsize',16)
axis square
box off
savefig(fullfile(figures_directory,'PV clustering.fig'))
saveas(gcf,fullfile(figures_directory,'PV clustering'),'png')

idx=kmeans(reduced_PV_correlations(:,end-2:end),max(reshaped_maps_vector));

figure
subplot(1,2,1)
for n=1:num_sessions*num_trials
    this_session=ceil(n/num_trials);
    if reshaped_maps_vector(n)>0
        plot(-reduced_PV_correlations(n,num_sessions*num_trials-1),-reduced_PV_correlations(n,num_sessions*num_trials),'.','markersize',20,'color',color_vector(reshaped_maps_vector(n),:))
        hold on
    end
end
xlabel('Component 1')
ylabel('Component 2')
set(gca,'fontsize',16)
axis square
box off
subplot(1,2,2)
for n=1:num_sessions*num_trials
    this_session=ceil(n/num_trials);
    if reshaped_maps_vector(n)>0
        plot(-reduced_PV_correlations(n,num_sessions*num_trials-1),-reduced_PV_correlations(n,num_sessions*num_trials),'.','markersize',20,'color',color_vector(idx(n),:))
        hold on
    end
end
xlabel('Component 1')
ylabel('Component 2')
set(gca,'fontsize',16)
axis square
box off

clustering_accuracy=zeros(1,factorial(max(reshaped_maps_vector)));
all_permutations=perms(1:max(reshaped_maps_vector));
for n=1:factorial(max(reshaped_maps_vector))
    this_permuted_clustering=zeros(size(idx));
    for k=1:max(reshaped_maps_vector)
        this_permuted_clustering(idx==k)=all_permutations(n,k);
    end
    clustering_accuracy(n)=100*sum(this_permuted_clustering(reshaped_maps_vector>0)==reshaped_maps_vector(reshaped_maps_vector>0))/sum(reshaped_maps_vector>0);
end
validated_clustering_accuracy=max(clustering_accuracy);

% long-term_dynamics:
within_maps_PV_dynamics=cell(1,num_sessions);
within_maps_separated_PV_dynamics=cell(1,num_sessions);
across_maps_separated_PV_dynamics=cell(1,num_sessions);
shuffle_separated_PV_dynamics=cell(1,num_sessions);
across_maps_flipped_PV_dynamics=cell(1,num_sessions);
across_maps_PV_dynamics=cell(1,num_sessions);
all_maps_PV_dynamics=cell(1,num_sessions);
within_maps_positional_shift_dynamics=cell(1,num_sessions);
across_maps_positional_shift_dynamics=cell(1,num_sessions);
all_maps_positional_shift_dynamics=cell(1,num_sessions);
within_maps_place_rate_dynamics=cell(1,num_sessions);
across_maps_place_rate_dynamics=cell(1,num_sessions);
within_maps_non_place_rate_dynamics=cell(1,num_sessions);
across_maps_non_place_rate_dynamics=cell(1,num_sessions);
all_rate_correlations_place_between_sessions=[];
all_rate_correlations_place_within_sessions=[];
all_rate_correlations_non_place_between_sessions=[];
all_rate_correlations_non_place_within_sessions=[];
same_map_count=0;
different_maps_count=0;
for n=1:num_sessions*num_trials
    first_session=ceil(n/num_trials);
    first_map=reshaped_maps_vector(n);
    for k=n+1:num_sessions*num_trials
        second_session=ceil(k/num_trials);
        second_map=reshaped_maps_vector(k);
        if first_session==second_session
            all_rate_correlations_place_within_sessions=[all_rate_correlations_place_within_sessions,rate_correlations_place_cells(n,k)];
            all_rate_correlations_non_place_within_sessions=[all_rate_correlations_non_place_within_sessions,rate_correlations_non_place_cells(n,k)];
        else
            all_rate_correlations_place_between_sessions=[all_rate_correlations_place_between_sessions,rate_correlations_place_cells(n,k)];
            all_rate_correlations_non_place_between_sessions=[all_rate_correlations_non_place_between_sessions,rate_correlations_non_place_cells(n,k)];
        end
        all_maps_PV_dynamics{second_session-first_session+1}=[all_maps_PV_dynamics{second_session-first_session+1},PV_correlations(n,k)];
        all_maps_positional_shift_dynamics{second_session-first_session+1}=[all_maps_positional_shift_dynamics{second_session-first_session+1},positional_shifts{n,k}];
        if first_map==second_map
            same_map_count=same_map_count+1;
            within_maps_PV_dynamics{second_session-first_session+1}=[within_maps_PV_dynamics{second_session-first_session+1},PV_correlations(n,k)];
            within_maps_positional_shift_dynamics{second_session-first_session+1}=[within_maps_positional_shift_dynamics{second_session-first_session+1},positional_shifts{n,k}];
            within_maps_place_rate_dynamics{second_session-first_session+1}=[within_maps_place_rate_dynamics{second_session-first_session+1},rate_correlations_place_cells(n,k)];
            within_maps_non_place_rate_dynamics{second_session-first_session+1}=[within_maps_non_place_rate_dynamics{second_session-first_session+1},rate_correlations_non_place_cells(n,k)];
            within_maps_separated_PV_dynamics{second_session-first_session+1}=[within_maps_separated_PV_dynamics{second_session-first_session+1},PV_correlations_separated(n,k)];
            shuffle_separated_PV_dynamics{second_session-first_session+1}=[shuffle_separated_PV_dynamics{second_session-first_session+1},shuffle_PV_correlations_separated(n,k)];
        else
            different_maps_count=different_maps_count+1;
            across_maps_PV_dynamics{second_session-first_session+1}=[across_maps_PV_dynamics{second_session-first_session+1},PV_correlations(n,k)];
            across_maps_positional_shift_dynamics{second_session-first_session+1}=[across_maps_positional_shift_dynamics{second_session-first_session+1},positional_shifts{n,k}];
            across_maps_place_rate_dynamics{second_session-first_session+1}=[across_maps_place_rate_dynamics{second_session-first_session+1},rate_correlations_place_cells(n,k)];
            across_maps_non_place_rate_dynamics{second_session-first_session+1}=[across_maps_non_place_rate_dynamics{second_session-first_session+1},rate_correlations_non_place_cells(n,k)];
            across_maps_separated_PV_dynamics{second_session-first_session+1}=[across_maps_separated_PV_dynamics{second_session-first_session+1},PV_correlations_separated(n,k)];
            across_maps_flipped_PV_dynamics{second_session-first_session+1}=[across_maps_flipped_PV_dynamics{second_session-first_session+1},flipped_PV_correlations(n,k)];
        end
    end
end
mean_within_maps_separated_PV_dynamics=zeros(1,num_sessions);
std_within_maps_separated_PV_dynamics=zeros(1,num_sessions);
mean_across_maps_separated_PV_dynamics=zeros(1,num_sessions);
std_across_maps_separated_PV_dynamics=zeros(1,num_sessions);
mean_across_maps_flipped_PV_dynamics=zeros(1,num_sessions);
std_across_maps_flipped_PV_dynamics=zeros(1,num_sessions);
mean_within_maps_PV_dynamics=zeros(1,num_sessions);
mean_across_maps_PV_dynamics=zeros(1,num_sessions);
mean_all_maps_PV_dynamics=zeros(1,num_sessions);
std_within_maps_PV_dynamics=zeros(1,num_sessions);
std_across_maps_PV_dynamics=zeros(1,num_sessions);
std_all_maps_PV_dynamics=zeros(1,num_sessions);
mean_within_maps_place_rate_dynamics=zeros(1,num_sessions);
mean_across_maps_place_rate_dynamics=zeros(1,num_sessions);
std_within_maps_place_rate_dynamics=zeros(1,num_sessions);
std_across_maps_place_rate_dynamics=zeros(1,num_sessions);
mean_within_maps_non_place_rate_dynamics=zeros(1,num_sessions);
mean_across_maps_non_place_rate_dynamics=zeros(1,num_sessions);
std_within_maps_non_place_rate_dynamics=zeros(1,num_sessions);
std_across_maps_non_place_rate_dynamics=zeros(1,num_sessions);
for n=1:num_sessions
    mean_within_maps_separated_PV_dynamics(n)=mean(within_maps_separated_PV_dynamics{n},'omitnan');
    std_within_maps_separated_PV_dynamics(n)=std(within_maps_separated_PV_dynamics{n},'omitnan');
    mean_across_maps_separated_PV_dynamics(n)=mean(across_maps_separated_PV_dynamics{n},'omitnan');
    std_across_maps_separated_PV_dynamics(n)=std(across_maps_separated_PV_dynamics{n},'omitnan');
    mean_across_maps_flipped_PV_dynamics(n)=mean(across_maps_flipped_PV_dynamics{n},'omitnan');
    std_across_maps_flipped_PV_dynamics(n)=std(across_maps_flipped_PV_dynamics{n},'omitnan');
    mean_within_maps_PV_dynamics(n)=mean(within_maps_PV_dynamics{n},'omitnan');
    std_within_maps_PV_dynamics(n)=std(within_maps_PV_dynamics{n},'omitnan');
    mean_across_maps_PV_dynamics(n)=mean(across_maps_PV_dynamics{n},'omitnan');
    std_across_maps_PV_dynamics(n)=std(across_maps_PV_dynamics{n},'omitnan');
    mean_all_maps_PV_dynamics(n)=mean(all_maps_PV_dynamics{n},'omitnan');
    std_all_maps_PV_dynamics(n)=std(all_maps_PV_dynamics{n},'omitnan');
    mean_within_maps_place_rate_dynamics(n)=mean(within_maps_place_rate_dynamics{n},'omitnan');
    std_within_maps_place_rate_dynamics(n)=std(within_maps_place_rate_dynamics{n},'omitnan');
    mean_across_maps_place_rate_dynamics(n)=mean(across_maps_place_rate_dynamics{n},'omitnan');
    std_across_maps_place_rate_dynamics(n)=std(across_maps_place_rate_dynamics{n},'omitnan');
    mean_within_maps_non_place_rate_dynamics(n)=mean(within_maps_non_place_rate_dynamics{n},'omitnan');
    std_within_maps_non_place_rate_dynamics(n)=std(within_maps_non_place_rate_dynamics{n},'omitnan');
    mean_across_maps_non_place_rate_dynamics(n)=mean(across_maps_non_place_rate_dynamics{n},'omitnan');
    std_across_maps_non_place_rate_dynamics(n)=std(across_maps_non_place_rate_dynamics{n},'omitnan');
end

% reorder PV correlations:
ordered_PV_correlations=zeros(sum(map_number_per_trial(:)>0),sum(map_number_per_trial(:)>0));
map_number_per_trial_transopsed=map_number_per_trial';
maps_indexes=map_number_per_trial_transopsed(:);
first_map=find(maps_indexes==1);
second_map=find(maps_indexes==2);
if number_of_maps>2
    third_map=find(maps_indexes==3);
end
if number_of_maps>3
    fourth_map=find(maps_indexes==4);
end
new_maps_order=[first_map;second_map];
if number_of_maps>2
    new_maps_order=[new_maps_order;third_map];
end
if number_of_maps>3
    new_maps_order=[new_maps_order;fourth_map];
end
for n=1:sum(map_number_per_trial(:)>0)
    for k=1:sum(map_number_per_trial(:)>0)
        ordered_PV_correlations(n,k)=PV_correlations(new_maps_order(n),new_maps_order(k));
    end
end

% plotting PV long-term dynamics:
figure
errorbar(2*(0:num_sessions-1),mean_within_maps_PV_dynamics,std_within_maps_PV_dynamics,'--*b')
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_PV_dynamics,std_across_maps_PV_dynamics,'--*r')
hold on
errorbar(2*(0:num_sessions-1),mean_all_maps_PV_dynamics,std_all_maps_PV_dynamics,'--*k')
xlim([0 2*(num_sessions-1)])
xlabel('Elapsed time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
legend('Within map','Across maps','All maps')
legend('boxoff')

figure
errorbar(2*(0:num_sessions-1),mean_within_maps_separated_PV_dynamics,std_within_maps_separated_PV_dynamics,'--*b')
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_separated_PV_dynamics,std_across_maps_separated_PV_dynamics,'--*r')
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_flipped_PV_dynamics,std_across_maps_flipped_PV_dynamics,'--*g')
xlim([0 2*(num_sessions-1)])
xlabel('Elapsed time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
legend('Within map','Across maps (same)','Across maps (flipped)')
legend('boxoff')

% plotting rate dynamics for place and non-place cells:
figure('units','normalized','outerposition',[0.3 0.3 0.4 0.4])
axes('position',[0.1 0.2 0.38 0.7])
errorbar(2*(0:num_sessions-1),mean_within_maps_place_rate_dynamics,std_within_maps_place_rate_dynamics,'--*b','linewidth',2)
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_place_rate_dynamics,std_across_maps_place_rate_dynamics,'--*r','linewidth',2)
xlim([0 2*(num_sessions-1)])
xlabel('Elapsed time (days)')
ylabel('Rate correlation')
title('Place cells','fontweight','normal')
set(gca,'fontsize',16)
box off
axes('position',[0.61 0.2 0.38 0.7])
errorbar(2*(0:num_sessions-1),mean_within_maps_non_place_rate_dynamics,std_within_maps_non_place_rate_dynamics,'--*b','linewidth',2)
hold on
errorbar(2*(0:num_sessions-1),mean_across_maps_non_place_rate_dynamics,std_across_maps_non_place_rate_dynamics,'--*r','linewidth',2)
xlim([0 2*(num_sessions-1)])
xlabel('Elapsed time (days)')
ylabel('Rate correlation')
title('Non-place cells','fontweight','normal')
set(gca,'fontsize',16)
box off

figure('units','normalized','outerposition',[0.3 0.3 0.4 0.4])
axes('position',[0.1 0.2 0.38 0.7])
scatter(all_rate_correlations_place_between_sessions,all_rate_correlations_non_place_between_sessions)
hold on
plot([0 1],[0 1],'--','color','k')
xlim([0 1])
ylim([0 1])
xlabel('Place cells')
ylabel('Non-place cells')
set(gca,'fontsize',16)
axes('position',[0.61 0.2 0.38 0.7])
scatter(all_rate_correlations_place_within_sessions,all_rate_correlations_non_place_within_sessions)
hold on
plot([0 1],[0 1],'--','color','k')
xlim([0 1])
ylim([0 1])
xlabel('Place cells')
ylabel('Non-place cells')
set(gca,'fontsize',16)

% plotting rate correlations for place and non-place cells:
figure('units','normalized','outerposition',[0.15 0.2 0.65 0.6])
axes('position',[0.035 0.12 0.39 0.8])
imagesc(rate_correlations_place_cells)
axis square
xlabel('Session number')
ylabel('Session number')
title('Place cells','fontweight','normal')
min_correlation=min(min(rate_correlations_place_cells));
caxis([round(min_correlation*100)/100 1])
set(gca,'xtick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'xticklabel',1:num_sessions)
set(gca,'ytick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'yticklabel',1:num_sessions)
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:num_sessions-1
    hold on
    plot(0.5+[0 num_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','w')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 num_sessions*num_trials],'linewidth',2,'color','w')
end
axes('position',[0.435 0.125 0.0145 0.79])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)

axes('position',[0.545 0.12 0.39 0.8])
imagesc(rate_correlations_non_place_cells)
axis square
xlabel('Session number')
ylabel('Session number')
title('Non-place cells','fontweight','normal')
min_correlation=min(min(rate_correlations_non_place_cells));
caxis([round(min_correlation*100)/100 1])
set(gca,'xtick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'xticklabel',1:num_sessions)
set(gca,'ytick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'yticklabel',1:num_sessions)
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:num_sessions-1
    hold on
    plot(0.5+[0 num_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','w')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 num_sessions*num_trials],'linewidth',2,'color','w')
end
axes('position',[0.945 0.125 0.0145 0.79])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)

% plotting rate map correlations:
figure('units','normalized','outerposition',[0.35 0.2 0.3 0.6])
axes('position',[0.1 0.1 0.75 0.75])
imagesc(rate_map_correlations)
xlabel('Session number')
ylabel('Session number')
min_correlation=min(min(rate_map_correlations));
caxis([round(min_correlation*100)/100 1])
set(gca,'xtick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'xticklabel',1:num_sessions)
set(gca,'ytick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'yticklabel',1:num_sessions)
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:num_sessions-1
    hold on
    plot(0.5+[0 num_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','w')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 num_sessions*num_trials],'linewidth',2,'color','w')
end
axes('position',[0.87 0.13 0.03 0.69])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate map correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)

figure
hist(rate_map_correlations(:),30)
xlabel('Rate map correlation')
ylabel('Count')
set(gca,'fontsize',16)

figure('units','normalized','outerposition',[0.35 0.2 0.3 0.6])
axes('position',[0.1 0.1 0.75 0.75])
imagesc(PV_correlations)
xlabel('Session number')
ylabel('Session number')
min_correlation=min(min(PV_correlations));
caxis([round(min_correlation*100)/100 1])
set(gca,'xtick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'xticklabel',1:num_sessions)
set(gca,'ytick',ceil(num_trials/2):num_trials:num_sessions*num_trials)
set(gca,'yticklabel',1:num_sessions)
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:num_sessions-1
    hold on
    plot(0.5+[0 num_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','w')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 num_sessions*num_trials],'linewidth',2,'color','w')
end
axes('position',[0.87 0.13 0.03 0.69])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'PV correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)
savefig(fullfile(figures_directory,'PV correlations.fig'))
saveas(gcf,fullfile(figures_directory,'PV correlations'),'png')

figure
imagesc(ordered_PV_correlations)
min_correlation=min(min(PV_correlations));
caxis([round(min_correlation*100)/100 1])
axis square
colormap('jet')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'fontsize',16)
savefig(fullfile(figures_directory,'Ordered PV correlations.fig'))
saveas(gcf,fullfile(figures_directory,'Ordered PV correlations'),'png')

figure
hist(PV_correlations(:),30)
xlabel('PV correlation')
ylabel('Count')
set(gca,'fontsize',16)

figure('units','normalized','outerposition',[0.35 0.2 0.3 0.6])
axes('position',[0.1 0.1 0.75 0.75])
imagesc(PV_correlations_bin_resolution)
colormap('jet')
set(gca,'fontsize',16)
savefig(fullfile(figures_directory,'PV correlation bin resolution.fig'))
saveas(gcf,fullfile(figures_directory,'PV correlation bin resolution'),'png')

selected_sessions=1:2;
num_selected_sessions=length(selected_sessions);
figure('units','normalized','outerposition',[0.35 0.2 0.3 0.6])
axes('position',[0.1 0.1 0.75 0.75])
imagesc(PV_correlations_bin_resolution((selected_sessions(1)-1)*num_trials*length(bins_to_use)+1:selected_sessions(end)*num_trials*length(bins_to_use),(selected_sessions(1)-1)*num_trials*length(bins_to_use)+1:selected_sessions(end)*num_trials*length(bins_to_use)))
set(gca,'xtick',[])
set(gca,'ytick',[])
min_correlation=round(min(min(PV_correlations_bin_resolution))*100)/100;
caxis([min_correlation 1])
colormap('jet')
axis square
for n=0:num_selected_sessions
    for k=1:num_trials-1
        hold on
        plot(0.5+[0 num_selected_sessions*num_trials*length(bins_to_use)+k*length(bins_to_use)],0.5+[n*num_trials*length(bins_to_use)+k*length(bins_to_use) n*num_trials*length(bins_to_use)+k*length(bins_to_use)],'linewidth',2,'color','k')
        hold on
        plot(0.5+[n*num_trials*length(bins_to_use)+k*length(bins_to_use) n*num_trials*length(bins_to_use)+k*length(bins_to_use)],0.5+[0 num_selected_sessions*num_trials*length(bins_to_use)+k*length(bins_to_use)],'linewidth',2,'color','k')
    end
end
for n=1:num_selected_sessions-1
    hold on
    plot(0.5+[0 num_selected_sessions*num_trials*length(bins_to_use)],0.5+[n*num_trials*length(bins_to_use) n*num_trials*length(bins_to_use)],'linewidth',3,'color','w')
    hold on
    plot(0.5+[n*num_trials*length(bins_to_use) n*num_trials*length(bins_to_use)],0.5+[0 num_selected_sessions*num_trials*length(bins_to_use)],'linewidth',3,'color','w')
end
axes('position',[0.87 0.13 0.03 0.69])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'PV correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(min_correlation),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)
savefig(fullfile(figures_directory,'PV correlation chosen sessions.fig'))
saveas(gcf,fullfile(figures_directory,'PV correlation chosen sessions'),'png')

figure
imagesc(PV_correlations_bin_resolution_separated)
colormap('jet')
savefig(fullfile(figures_directory,'PV correlation bin resolution - separated.fig'))
saveas(gcf,fullfile(figures_directory,'PV correlation bin resolution - separated'),'png')

%% Calculating maps and correlations within a chosen session:
chosen_session=2;
chosen_trials=[1 5];
maps_vector=map_number_per_trial(chosen_session,:);

% calculating rate map and PV correlations:
rate_map_correlations_chosen_session=ones(num_trials,num_trials);
PV_correlations_chosen_session=ones(num_trials,num_trials);
PV_correlations_bin_resolution_chosen_session=zeros(num_trials*length(bins_to_use),num_trials*length(bins_to_use));
all_rate_map_correlations_same=[];
all_rate_map_correlations_diff=[];
shuffle_rate_map_correlations=[];
shuffle_rate_map_correlations_ver2=[];
shuffle_rate_map_correlations_ver3=[];
all_rate_map_correlations_same_flipped=[];
all_rate_map_correlations_diff_flipped=[];
all_PV_correlations_same=[];
all_PV_correlations_diff=[];
shuffle_PV_correlations=[];
all_PV_correlations_same_flipped=[];
all_PV_correlations_diff_flipped=[];

for s=1:num_sessions
    for n=1:num_trials
        for k=n:num_trials
            this_significant=significant_sessions{s};
            all_maps_trial_one=squeeze(smoothed_rate_maps_trials_both_sides(s,n,this_significant,:));
            all_maps_trial_two=squeeze(smoothed_rate_maps_trials_both_sides(s,k,this_significant,:));
            
            [~,random_indexes]=sort(rand(1,size(all_maps_trial_one,1)));
            shuffle_all_maps_trial_one=all_maps_trial_one(random_indexes,:);
            
            [~,random_indexes]=sort(rand(1,size(all_maps_trial_one,2)));
            shuffle_all_maps_trial_one_ver2=all_maps_trial_one(:,random_indexes);
                        
            [~,random_indexes]=sort(rand(size(all_maps_trial_one,2),size(all_maps_trial_one,1)));
            shuffle_all_maps_trial_one_ver3=all_maps_trial_one(random_indexes');
          
            all_rate_map_correlations_temp=corr(all_maps_trial_one',all_maps_trial_two');
            shuffle_all_rate_map_correlations_temp=corr(shuffle_all_maps_trial_one',all_maps_trial_two');
            shuffle_all_rate_map_correlations_temp_ver2=corr(shuffle_all_maps_trial_one_ver2',all_maps_trial_two');
            shuffle_all_rate_map_correlations_temp_ver3=corr(shuffle_all_maps_trial_one_ver3',all_maps_trial_two');
            all_rate_map_correlations_flipped=corr(fliplr(all_maps_trial_one)',all_maps_trial_two');
            rate_map_correlations_temp=zeros(1,length(this_significant));
            shuffle_rate_map_correlations_temp=zeros(1,length(this_significant));
            shuffle_rate_map_correlations_temp_ver2=zeros(1,length(this_significant));
            shuffle_rate_map_correlations_temp_ver3=zeros(1,length(this_significant));
            rate_map_correlations_flipped_temp=zeros(1,length(this_significant));
            for m=1:length(this_significant)
                rate_map_correlations_temp(m)=all_rate_map_correlations_temp(m,m);
                shuffle_rate_map_correlations_temp(m)=shuffle_all_rate_map_correlations_temp(m,m);
                shuffle_rate_map_correlations_temp_ver2(m)=shuffle_all_rate_map_correlations_temp_ver2(m,m);
                shuffle_rate_map_correlations_temp_ver3(m)=shuffle_all_rate_map_correlations_temp_ver3(m,m);
                rate_map_correlations_flipped_temp(m)=all_rate_map_correlations_flipped(m,m);
            end
            all_PV_correlations_temp=corr(all_maps_trial_one,all_maps_trial_two);
            shuffle_all_PV_correlations_temp=corr(shuffle_all_maps_trial_one,all_maps_trial_two);
            all_PV_correlations_flipped_temp=corr(fliplr(all_maps_trial_one),all_maps_trial_two);
            PV_correlations_temp=zeros(1,length(bins_to_use));
            shuffle_PV_correlations_temp=zeros(1,length(bins_to_use));
            PV_correlations_flipped_temp=zeros(1,length(bins_to_use));
            if k~=n
                if map_number_per_trial(s,n)==map_number_per_trial(s,k)
                    all_rate_map_correlations_same=[all_rate_map_correlations_same , rate_map_correlations_temp(~isnan(rate_map_correlations_temp))];
                    shuffle_rate_map_correlations=[shuffle_rate_map_correlations , shuffle_rate_map_correlations_temp(~isnan(shuffle_rate_map_correlations_temp))];
                    shuffle_rate_map_correlations_ver2=[shuffle_rate_map_correlations_ver2 , shuffle_rate_map_correlations_temp(~isnan(shuffle_rate_map_correlations_temp_ver2))];
                    shuffle_rate_map_correlations_ver3=[shuffle_rate_map_correlations_ver3 , shuffle_rate_map_correlations_temp(~isnan(shuffle_rate_map_correlations_temp_ver3))];
                    all_rate_map_correlations_same_flipped=[all_rate_map_correlations_same_flipped , rate_map_correlations_flipped_temp(~isnan(rate_map_correlations_flipped_temp))];
                else
                    all_rate_map_correlations_diff=[all_rate_map_correlations_diff , rate_map_correlations_temp(~isnan(rate_map_correlations_temp))];
                    all_rate_map_correlations_diff_flipped=[all_rate_map_correlations_diff_flipped , rate_map_correlations_flipped_temp(~isnan(rate_map_correlations_flipped_temp))];
                end
                for m=1:length(bins_to_use)
                    PV_correlations_temp(m)=all_PV_correlations_temp(m,m);
                    PV_correlations_flipped_temp(m)=all_PV_correlations_flipped_temp(m,m);
                    shuffle_PV_correlations_temp(m)=shuffle_all_PV_correlations_temp(m,m);
                end
                if map_number_per_trial(s,n)==map_number_per_trial(s,k)
                    all_PV_correlations_same=[all_PV_correlations_same , PV_correlations_temp(~isnan(PV_correlations_temp))];
                    all_PV_correlations_same_flipped=[all_PV_correlations_same_flipped , PV_correlations_flipped_temp(~isnan(PV_correlations_flipped_temp))];
                    shuffle_PV_correlations=[shuffle_PV_correlations, shuffle_PV_correlations_temp(~isnan(shuffle_PV_correlations_temp))];
                else
                    all_PV_correlations_diff=[all_PV_correlations_diff , PV_correlations_temp(~isnan(PV_correlations_temp))];
                    all_PV_correlations_diff_flipped=[all_PV_correlations_diff_flipped , PV_correlations_flipped_temp(~isnan(PV_correlations_flipped_temp))];
                end
                if chosen_session==s
                    rate_map_correlations_chosen_session(n,k)=mean(rate_map_correlations_temp,'omitnan');
                    rate_map_correlations_chosen_session(k,n)=rate_map_correlations_chosen_session(n,k);
                    PV_correlations_chosen_session(n,k)=mean(PV_correlations_temp,'omitnan');
                    PV_correlations_chosen_session(k,n)=PV_correlations_chosen_session(n,k);
                end
            end
            if chosen_session==s
                PV_correlations_bin_resolution_chosen_session((n-1)*length(bins_to_use)+1:n*length(bins_to_use),(k-1)*length(bins_to_use)+1:k*length(bins_to_use))=all_PV_correlations_temp;
                PV_correlations_bin_resolution_chosen_session((k-1)*length(bins_to_use)+1:k*length(bins_to_use),(n-1)*length(bins_to_use)+1:n*length(bins_to_use))=all_PV_correlations_temp';
            end
        end
    end
end

% ordering the place fields:
significant_chosen_session=significant_sessions{chosen_session};
significant_chosen_session(significant_chosen_session>number_of_neurons)=significant_chosen_session(significant_chosen_session>number_of_neurons)-number_of_neurons;
significant_chosen_session=unique(significant_chosen_session);
[map_1_values,ordered_locations_map_1]=sort(preferred_locations_averaged_both_sides(chosen_session,chosen_trials(1),significant_chosen_session));
[map_2_values,ordered_locations_map_2]=sort(preferred_locations_averaged_both_sides(chosen_session,chosen_trials(2),significant_chosen_session));

ordered_normalized_map_1=zeros(num_trials,length(significant_chosen_session),length(bins_to_use));
corresponding_order_map_2=zeros(num_trials,length(significant_chosen_session),length(bins_to_use));
for m=1:num_trials
    for n=1:length(significant_chosen_session)
         corresponding_order_map_2(m,n,:)=normalized_rate_maps_averaged_both_sides(chosen_session,m,significant_chosen_session(ordered_locations_map_1(n)),:);
        ordered_normalized_map_1(m,n,:)=normalized_rate_maps_averaged_both_sides(chosen_session,m,significant_chosen_session(ordered_locations_map_1(n)),:);
    end
end
ordered_normalized_map_2=zeros(num_trials,length(significant_chosen_session),length(bins_to_use));
corresponding_order_map_1=zeros(num_trials,length(significant_chosen_session),length(bins_to_use));
for m=1:num_trials
    for n=1:length(significant_chosen_session)
        corresponding_order_map_1(m,n,:)=normalized_rate_maps_averaged_both_sides(chosen_session,m,significant_chosen_session(ordered_locations_map_2(n)),:);
        ordered_normalized_map_2(m,n,:)=normalized_rate_maps_averaged_both_sides(chosen_session,m,significant_chosen_session(ordered_locations_map_2(n)),:);
    end
end

% calculating the PV correlations between consecutive time windows:
smoothing_kernel=gausswin(1+round(length(bins_to_use)/4));
time_window=60;
correlation_resolution=10;
significant_cells=significant_sessions{chosen_session};
significant_cells(significant_cells>number_of_neurons)=significant_cells(significant_cells>number_of_neurons)-number_of_neurons;
significant_cells=unique(significant_cells);
this_comparison_number_of_neurons=length(significant_cells);
auto_correlation_vector=zeros(1,floor((size(concatanated_binary_events{chosen_session},1)-fs*2*time_window)/(fs*correlation_resolution)));
t_auto_correlation=time_window:correlation_resolution:time_window+correlation_resolution*(length(auto_correlation_vector)-1);
for n=1:length(auto_correlation_vector)
    before_indexes=(t_auto_correlation(n)-time_window)*fs+1:t_auto_correlation(n)*fs;
    after_indexes=t_auto_correlation(n)*fs+1:(t_auto_correlation(n)+time_window)*fs;
    
    first_trial_index=find(before_indexes(100)<cumsum(number_of_frames_per_trial(chosen_session,:)),1,'first');
    second_trial_index=find(after_indexes(end-100)<cumsum(number_of_frames_per_trial(chosen_session,:)),1,'first');
    rate_maps_first_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
    rate_maps_second_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
    smoothed_rate_maps_first_window=zeros(size(rate_maps_first_window));
    smoothed_rate_maps_second_window=zeros(size(rate_maps_second_window));
    
    for k=1:length(bins_to_use)
        this_bins_right_index_before=find(round(concatanated_locations{chosen_session}(before_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(before_indexes)>velocity_thresh);
        this_bins_left_index_before=find(round(concatanated_locations{chosen_session}(before_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(before_indexes)<-velocity_thresh);
        this_bins_right_index_after=find(round(concatanated_locations{chosen_session}(after_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(after_indexes)>velocity_thresh);
        this_bins_left_index_after=find(round(concatanated_locations{chosen_session}(after_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(after_indexes)<-velocity_thresh);
        rate_maps_first_window(:,k,1)=mean(concatanated_binary_events{chosen_session}(before_indexes(1)-1+this_bins_right_index_before,significant_cells));
        rate_maps_first_window(:,k,2)=mean(concatanated_binary_events{chosen_session}(before_indexes(1)-1+this_bins_left_index_before,significant_cells));
        rate_maps_second_window(:,k,1)=mean(concatanated_binary_events{chosen_session}(after_indexes(1)-1+this_bins_right_index_after,significant_cells));
        rate_maps_second_window(:,k,2)=mean(concatanated_binary_events{chosen_session}(after_indexes(1)-1+this_bins_left_index_after,significant_cells));
    end
    
    for m=1:this_comparison_number_of_neurons
        smoothed_rate_maps_first_window(m,:,1)=conv(rate_maps_first_window(m,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        smoothed_rate_maps_first_window(m,:,2)=conv(rate_maps_first_window(m,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        smoothed_rate_maps_second_window(m,:,1)=conv(rate_maps_second_window(m,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        smoothed_rate_maps_second_window(m,:,2)=conv(rate_maps_second_window(m,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
    end
    
    smoothed_rate_maps_first_window_both_sides=[squeeze(smoothed_rate_maps_first_window(:,:,1)) ; squeeze(smoothed_rate_maps_first_window(:,:,2))];
    smoothed_rate_maps_second_window_both_sides=[squeeze(smoothed_rate_maps_second_window(:,:,1)) ; squeeze(smoothed_rate_maps_second_window(:,:,2))];
    
    all_PV_correlations_temp=corr(smoothed_rate_maps_first_window_both_sides,smoothed_rate_maps_second_window_both_sides);
    PV_correlations_temp=zeros(1,length(bins_to_use));
    for k=1:length(bins_to_use)
        PV_correlations_temp(k)=all_PV_correlations_temp(k,k);
    end
    auto_correlation_vector(n)=mean(PV_correlations_temp,'omitnan');
end

% calculating within trial correlations:
within_trials_PV_correlations=zeros(1,num_sessions*num_trials);
trial_count=0;
for n=1:num_sessions
    significant_cells=significant_sessions{n};
    significant_cells(significant_cells>number_of_neurons)=significant_cells(significant_cells>number_of_neurons)-number_of_neurons;
    significant_cells=unique(significant_cells);
    this_comparison_number_of_neurons=length(significant_cells);
    for m=1:num_trials
        trial_count=trial_count+1;
        this_trial_length=number_of_frames_per_trial(n,m);
        first_half_indexes=1:round(this_trial_length/2);
        second_half_indexes=round(this_trial_length/2)+1:this_trial_length;
        
        rate_maps_first_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
        rate_maps_second_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
        smoothed_rate_maps_first_window=zeros(size(rate_maps_first_window));
        smoothed_rate_maps_second_window=zeros(size(rate_maps_second_window));
        
        for k=1:length(bins_to_use)
            this_bins_right_index_before=find(round(all_locations{n,m}(first_half_indexes))==bins_to_use(k) & all_velocities{n,m}(first_half_indexes)>velocity_thresh);
            this_bins_left_index_before=find(round(all_locations{n,m}(first_half_indexes))==bins_to_use(k) & all_velocities{n,m}(first_half_indexes)<-velocity_thresh);
            this_bins_right_index_after=find(round(all_locations{n,m}(second_half_indexes))==bins_to_use(k) & all_velocities{n,m}(second_half_indexes)>velocity_thresh);
            this_bins_left_index_after=find(round(all_locations{n,m}(second_half_indexes))==bins_to_use(k) & all_velocities{n,m}(second_half_indexes)<-velocity_thresh);
            rate_maps_first_window(:,k,1)=mean(all_events{n,m}(first_half_indexes(1)-1+this_bins_right_index_before,significant_cells));
            rate_maps_first_window(:,k,2)=mean(all_events{n,m}(first_half_indexes(1)-1+this_bins_left_index_before,significant_cells));
            rate_maps_second_window(:,k,1)=mean(all_events{n,m}(second_half_indexes(1)-1+this_bins_right_index_after,significant_cells));
            rate_maps_second_window(:,k,2)=mean(all_events{n,m}(second_half_indexes(1)-1+this_bins_left_index_after,significant_cells));
        end
        
        for k=1:this_comparison_number_of_neurons
            smoothed_rate_maps_first_window(k,:,1)=conv(rate_maps_first_window(k,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_first_window(k,:,2)=conv(rate_maps_first_window(k,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_second_window(k,:,1)=conv(rate_maps_second_window(k,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_second_window(k,:,2)=conv(rate_maps_second_window(k,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        end
        
        smoothed_rate_maps_first_window_both_sides=[squeeze(smoothed_rate_maps_first_window(:,:,1)) ; squeeze(smoothed_rate_maps_first_window(:,:,2))];
        smoothed_rate_maps_second_window_both_sides=[squeeze(smoothed_rate_maps_second_window(:,:,1)) ; squeeze(smoothed_rate_maps_second_window(:,:,2))];
        
        all_PV_correlations_temp=corr(smoothed_rate_maps_first_window_both_sides,smoothed_rate_maps_second_window_both_sides);
        PV_correlations_temp=zeros(1,length(bins_to_use));
        for k=1:length(bins_to_use)
            PV_correlations_temp(k)=all_PV_correlations_temp(k,k);
        end
        within_trials_PV_correlations(trial_count)=mean(PV_correlations_temp,'omitnan');
    end
end

% plotting within trial correlations
figure('units','normalized','outerposition',[0.3 0.3 0.4 0.4])
axes('position',[0.1 0.2 0.38 0.75])
plot(t_auto_correlation,auto_correlation_vector,'linewidth',2)
xlabel('Time (sec)')
ylabel('PV correlation')
for n=1:num_trials-1
    hold on
    plot([sum(number_of_frames_per_trial(chosen_session,1:n))/fs sum(number_of_frames_per_trial(chosen_session,1:n))/fs],[-0.2 1],'--','linewidth',2,'color','r')
    text(90+(n-1)*180,0.75,['T' num2str(n)],'HorizontalAlignment','center','fontsize',16)
end
text(90+(4)*180,0.75,'T5','HorizontalAlignment','center','fontsize',16)
set(gca,'fontsize',16)
ylim([-0.1 1])
xlim([0 900])
box off

axes('position',[0.61 0.2 0.38 0.75])
x_vec=0.025:0.05:0.975;
[within_trial_PV_distribution,~]=hist(within_trials_PV_correlations,x_vec);
within_trial_PV_distribution=within_trial_PV_distribution./sum(within_trial_PV_distribution);
bar(x_vec,within_trial_PV_distribution,1)
xlim([0 1])
xlabel('PV correlation')
ylabel('Frequency')
set(gca,'fontsize',16)
box off

% plotting the place field maps:
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
for m=1:num_trials
    axes('position',[(m-1)*0.185+0.07 0.6 0.15 0.35])
    imagesc(squeeze(ordered_normalized_map_1(m,~isnan(map_1_values),:)))
    color_map=colormap('jet');
    if m==1
        xlabel('Position (cm)')
        ylabel('Cell number')
        set(gca,'xtick',linspace(1,24,7))
        set(gca,'xticklabel',linspace(0,track_length,7))
    else
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    title(['Trial ' num2str(m)],'fontweight','normal')
    set(gca,'fontsize',12)
    if m==num_trials
        axes('position',[0.97 0.6 0.01 0.35])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        xlim([0 1])
        ylim([0 1])
        for n=1:size(color_map,1)
            hold on
            p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
            set(p,'FaceAlpha',1,'EdgeColor','none');
        end
        plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
        text(2,0.5,'Normalized rate','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
        text(1.2,0,'0','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
        text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
    end
    
    axes('position',[(m-1)*0.185+0.07 0.12 0.15 0.35])
    imagesc(squeeze(ordered_normalized_map_2(m,~isnan(map_2_values),:)))
    colormap('jet')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'fontsize',12)
end

% plotting PV correlations:
figure('units','normalized','outerposition',[0.2 0.3 0.6 0.4])
axes('position',[0.05 0.2 0.25 0.7])
imagesc(PV_correlations_bin_resolution_chosen_session)
min_correlation=min(min(PV_correlations_bin_resolution_chosen_session));
caxis([round(min_correlation*100)/100 1])
xlabel('Trial number')
ylabel('Trial number')
for n=1:num_trials-1
    hold on
    plot(0.5+[0 num_trials*length(bins_to_use)],0.5+[n*length(bins_to_use) n*length(bins_to_use)],'linewidth',2,'color','w')
    hold on
    plot(0.5+[n*length(bins_to_use) n*length(bins_to_use)],0.5+[0 num_trials*n*length(bins_to_use)],'linewidth',2,'color','w')
end

set(gca,'xtick',length(bins_to_use)/2:length(bins_to_use):length(bins_to_use)*num_trials)
set(gca,'xticklabel',1:num_trials)
set(gca,'ytick',length(bins_to_use)/2:length(bins_to_use):length(bins_to_use)*num_trials)
set(gca,'yticklabel',1:num_trials)
title('PV correlations','fontweight','normal')
set(gca,'fontsize',16)
axis square
colormap('jet')
axes('position',[0.29 0.2 0.01 0.7])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)

axes('position',[0.39 0.2 0.25 0.7])
x_vec=round(min(all_PV_correlations_diff)*100)/100:0.05:1;
[n1,~]=hist(all_PV_correlations_same,x_vec);
[n2,~]=hist(all_PV_correlations_diff,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
bar(x_vec,n1,'b','barwidth',1);
hold on
bar(x_vec,n2,'r','barwidth',1);
xlim([-0.2 1])
xlabel('PV correlation')
ylabel('Frequency')
set(gca,'fontsize',16)
legend('Same map','Different maps')
legend('boxoff')
text(-0.3,-0.4,'Frequency','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
box off

axes('position',[0.73 0.2 0.25 0.7])
x_vec=round(min(all_rate_map_correlations_diff)*100)/100:0.05:1;
[n1,~]=hist(all_rate_map_correlations_same,x_vec);
[n2,~]=hist(all_rate_map_correlations_diff,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
bar(x_vec,n1,'b','barwidth',1);
hold on
bar(x_vec,n2,'r','barwidth',1);
xlim([-0.6 1])
xlabel('Rate map correlation')
ylabel('Frequency')
set(gca,'fontsize',16)
legend('Same map','Different maps')
legend('boxoff')
text(-0.3,-0.4,'Frequency','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
box off

% plotting PV correlations for same versus flipped orientations:
figure('units','normalized','outerposition',[0.2 0.3 0.6 0.4])
axes('position',[0.39 0.2 0.25 0.7])
x_vec=round(min(all_PV_correlations_diff_flipped)*100)/100:0.05:1;
[n1,~]=hist(all_PV_correlations_same_flipped,x_vec);
[n2,~]=hist(all_PV_correlations_diff_flipped,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
bar(x_vec,n1,'b','barwidth',1);
hold on
bar(x_vec,n2,'r','barwidth',1);
xlim([-0.2 1])
xlabel('PV correlation')
ylabel('Frequency')
set(gca,'fontsize',16)
legend('Same map flipped','Different maps flipped')
legend('boxoff')
text(-0.3,-0.4,'Frequency','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
box off

axes('position',[0.73 0.2 0.25 0.7])
x_vec=round(min(all_rate_map_correlations_diff_flipped)*100)/100:0.05:1;
[n1,~]=hist(all_rate_map_correlations_same_flipped,x_vec);
[n2,~]=hist(all_rate_map_correlations_diff_flipped,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
bar(x_vec,n1,'b','barwidth',1);
hold on
bar(x_vec,n2,'r','barwidth',1);
xlim([-0.6 1])
xlabel('Rate map correlation')
ylabel('Frequency')
set(gca,'fontsize',16)
legend('Same map flipped','Different maps flipped')
legend('boxoff')
text(-0.3,-0.4,'Frequency','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
box off

%%
% calculating the PV correlations between consecutive time windows:
smoothing_kernel=gausswin(1+round(length(bins_to_use)/4));
smoothing_kernel_2=gausswin(5);
time_window=60;
correlation_resolution=10;
switch_threshold=0.1;
end_trial_window_size=30;

for p=1:num_sessions
    significant_cells=significant_sessions{p};
    significant_cells(significant_cells>number_of_neurons)=significant_cells(significant_cells>number_of_neurons)-number_of_neurons;
    significant_cells=unique(significant_cells);
    this_comparison_number_of_neurons=length(significant_cells);
    auto_correlation_vector=zeros(1,floor((size(concatanated_binary_events{p},1)-fs*2*time_window)/(fs*correlation_resolution)));
    time_of_trials=frame_log{p}(relevant_trials,3:4)/fs;
    time_of_trials_plus_window=time_of_trials+repmat([-end_trial_window_size/2 end_trial_window_size/2],num_trials,1);
    end_trial_window_range=zeros(length(relevant_trials)-1,2);
    for l=1:num_trials-1
        end_trial_window_range(l,1)=time_of_trials_plus_window(l+1,1);
        end_trial_window_range(l,2)=time_of_trials_plus_window(l,2);
    end
    t_auto_correlation=time_window:correlation_resolution:time_window+correlation_resolution*(length(auto_correlation_vector)-1);
    end_trial_window_range=end_trial_window_range-time_of_trials(1,1);
    
    for n=1:length(auto_correlation_vector)
        before_indexes=(t_auto_correlation(n)-time_window)*fs+1:t_auto_correlation(n)*fs;
        after_indexes=t_auto_correlation(n)*fs+1:(t_auto_correlation(n)+time_window)*fs;
        
        first_trial_index=find(before_indexes(100)<cumsum(number_of_frames_per_trial(p,:)),1,'first');
        second_trial_index=find(after_indexes(end-100)<cumsum(number_of_frames_per_trial(p,:)),1,'first');
        rate_maps_first_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
        rate_maps_second_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
        smoothed_rate_maps_first_window=zeros(size(rate_maps_first_window));
        smoothed_rate_maps_second_window=zeros(size(rate_maps_second_window));
        
        for k=1:length(bins_to_use)
            this_bins_right_index_before=find(round(concatanated_locations{p}(before_indexes))==bins_to_use(k) & concatanated_velocities{p}(before_indexes)>velocity_thresh);
            this_bins_left_index_before=find(round(concatanated_locations{p}(before_indexes))==bins_to_use(k) & concatanated_velocities{p}(before_indexes)<-velocity_thresh);
            this_bins_right_index_after=find(round(concatanated_locations{p}(after_indexes))==bins_to_use(k) & concatanated_velocities{p}(after_indexes)>velocity_thresh);
            this_bins_left_index_after=find(round(concatanated_locations{p}(after_indexes))==bins_to_use(k) & concatanated_velocities{p}(after_indexes)<-velocity_thresh);
            rate_maps_first_window(:,k,1)=mean(concatanated_binary_events{p}(before_indexes(1)-1+this_bins_right_index_before,significant_cells));
            rate_maps_first_window(:,k,2)=mean(concatanated_binary_events{p}(before_indexes(1)-1+this_bins_left_index_before,significant_cells));
            rate_maps_second_window(:,k,1)=mean(concatanated_binary_events{p}(after_indexes(1)-1+this_bins_right_index_after,significant_cells));
            rate_maps_second_window(:,k,2)=mean(concatanated_binary_events{p}(after_indexes(1)-1+this_bins_left_index_after,significant_cells));
        end
        
        for m=1:this_comparison_number_of_neurons
            smoothed_rate_maps_first_window(m,:,1)=conv(rate_maps_first_window(m,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_first_window(m,:,2)=conv(rate_maps_first_window(m,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_second_window(m,:,1)=conv(rate_maps_second_window(m,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_second_window(m,:,2)=conv(rate_maps_second_window(m,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        end
        
        smoothed_rate_maps_first_window_both_sides=[squeeze(smoothed_rate_maps_first_window(:,:,1)) ; squeeze(smoothed_rate_maps_first_window(:,:,2))];
        smoothed_rate_maps_second_window_both_sides=[squeeze(smoothed_rate_maps_second_window(:,:,1)) ; squeeze(smoothed_rate_maps_second_window(:,:,2))];
        
        all_PV_correlations_temp=corr(smoothed_rate_maps_first_window_both_sides,smoothed_rate_maps_second_window_both_sides);
        PV_correlations_temp=zeros(1,length(bins_to_use));
        for k=1:length(bins_to_use)
            PV_correlations_temp(k)=all_PV_correlations_temp(k,k);
        end
        auto_correlation_vector(n)=mean(PV_correlations_temp,'omitnan');
    end
end

%% Place field relations across maps:
averaged_smoothed_rate_maps=zeros(number_of_maps,number_of_neurons,length(bins_to_use),2);
averaged_smoothed_rate_maps_both_sides=zeros(number_of_maps,2*number_of_neurons,length(bins_to_use));
concatenated_smoothed_rate_maps=zeros(number_of_maps,number_of_neurons,2*length(bins_to_use));
for n=1:number_of_maps
    indexes_in_this_map=find(reshaped_maps_vector==n);
    sessions_in_this_map=ceil(indexes_in_this_map/num_trials);
    trials_in_this_map=indexes_in_this_map-(sessions_in_this_map-1)*num_trials;
    temp_rate_maps=zeros(length(indexes_in_this_map),number_of_neurons,length(bins_to_use),2);
    for k=1:length(indexes_in_this_map)
        temp_rate_maps(k,:,:,:)=smoothed_rate_maps_trials(sessions_in_this_map(k),trials_in_this_map(k),:,:,:);
    end
    if length(trials_in_this_map)>1
        averaged_smoothed_rate_maps(n,:,:,:)=mean(temp_rate_maps);
    else
        averaged_smoothed_rate_maps(n,:,:,:)=squeeze(temp_rate_maps);
    end
    averaged_smoothed_rate_maps_both_sides(n,:,:)=[squeeze(averaged_smoothed_rate_maps(n,:,:,1)) ; squeeze(averaged_smoothed_rate_maps(n,:,:,2))];
    concatenated_smoothed_rate_maps(n,:,:)=[squeeze(averaged_smoothed_rate_maps(n,:,:,1)) , squeeze(averaged_smoothed_rate_maps(n,:,:,2))];
end

all_significant_right=all_significant(all_significant<=number_of_neurons);
all_significant_left=all_significant(all_significant>number_of_neurons);
all_significant_left=all_significant_left-number_of_neurons;
all_significant_both_directions=unique(union(all_significant_left,all_significant_right));

PV_correlations_per_map_bin_resolution=zeros(number_of_maps*length(bins_to_use)*2,number_of_maps*length(bins_to_use)*2);
PV_correlations_per_map_both_sides=zeros(number_of_maps*length(bins_to_use),number_of_maps*length(bins_to_use));
rate_map_correlations_across_maps=cell(sum(1:number_of_maps-1),2*length(bins_to_use));
rate_map_correlations_across_maps_flipped=cell(sum(1:number_of_maps-1),2*length(bins_to_use));
rate_map_correlations_per_map=cell(1,number_of_maps);
all_concatenated_rate_map_correlations=cell(1,sum(1:number_of_maps-1));
comparison_count=0;
all_cells_predicted_stability=cell(1,sum(1:number_of_maps-1));
all_cells_stability=cell(1,sum(1:number_of_maps-1));
all_cells_predicted_stability_shuffle=cell(1,sum(1:number_of_maps-1));
for n=1:number_of_maps
    all_maps_first_map_right=squeeze(averaged_smoothed_rate_maps(n,all_significant_both_directions,:,1));
    all_maps_first_map_left=squeeze(averaged_smoothed_rate_maps(n,all_significant_both_directions,:,2));
    all_maps_first_map=[squeeze(all_maps_first_map_right) , squeeze(all_maps_first_map_left)];
    all_neurons_maps_first_map_right=squeeze(averaged_smoothed_rate_maps(n,:,:,1));
    all_neurons_maps_first_map_left=squeeze(averaged_smoothed_rate_maps(n,:,:,2));
    all_neurons_maps_first_map=[squeeze(all_neurons_maps_first_map_right) , squeeze(all_neurons_maps_first_map_left)];
    all_concatenated_maps_first_map=squeeze(concatenated_smoothed_rate_maps(n,:,:));
    [max_activity_first_map,preferred_locations_first_map]=max(all_maps_first_map');
    preferred_locations_first_map(max_activity_first_map==0)=nan;
    all_pv_correlations_temp=corr(all_maps_first_map,all_maps_first_map);
    all_pv_correlations_both_sides_temp=corr(squeeze(averaged_smoothed_rate_maps_both_sides(n,:,:)),squeeze(averaged_smoothed_rate_maps_both_sides(n,:,:)));
    PV_correlations_per_map_bin_resolution((n-1)*2*length(bins_to_use)+1:n*2*length(bins_to_use),(n-1)*2*length(bins_to_use)+1:n*2*length(bins_to_use))=all_pv_correlations_temp;
    PV_correlations_per_map_both_sides((n-1)*length(bins_to_use)+1:n*length(bins_to_use),(n-1)*length(bins_to_use)+1:n*length(bins_to_use))=all_pv_correlations_both_sides_temp;
    
    rate_map_correlations_per_map{n}=nan(number_of_neurons,number_of_neurons);
    rate_map_correlations_per_map{n}=corr(all_neurons_maps_first_map',all_neurons_maps_first_map');
    
    for k=n+1:number_of_maps
        comparison_count=comparison_count+1;
        cell_count=0;
        all_maps_second_map_right=squeeze(averaged_smoothed_rate_maps(k,all_significant_both_directions,:,1));
        all_maps_second_map_left=squeeze(averaged_smoothed_rate_maps(k,all_significant_both_directions,:,2));
        all_maps_second_map=[squeeze(all_maps_second_map_right) , squeeze(all_maps_second_map_left)];
        all_concatenated_maps_second_map=squeeze(concatenated_smoothed_rate_maps(k,:,:));
        
        [max_activity_second_map,preferred_locations_second_map]=max(all_maps_second_map');
        preferred_locations_second_map(max_activity_second_map==0)=nan;
        cells_with_place_field=find(~isnan(preferred_locations_first_map) & ~isnan(preferred_locations_second_map));
        all_cells_predicted_stability{comparison_count}=zeros(1,length(cells_with_place_field));
        all_cells_stability{comparison_count}=zeros(1,length(cells_with_place_field));
        all_cells_predicted_stability_shuffle{comparison_count}=zeros(1,length(cells_with_place_field));
        cells_count=0;
        for m=1:length(cells_with_place_field)
            this_bin=preferred_locations_first_map(cells_with_place_field(m));
            other_cells_this_bin=find(preferred_locations_first_map(cells_with_place_field)==this_bin);
            other_cells_this_bin=setdiff(other_cells_this_bin,m);
            if ~isempty(other_cells_this_bin)
                cells_count=cells_count+1;
                all_cells_predicted_stability{comparison_count}(cells_count)=trace(corr(all_maps_first_map(cells_with_place_field(other_cells_this_bin),:)',fliplr(all_maps_second_map(cells_with_place_field(other_cells_this_bin),:))'))/length(other_cells_this_bin);
                [~,shuffle_cells_this_bin]=sort(rand(1,length(cells_with_place_field)));
                all_cells_predicted_stability_shuffle{comparison_count}(cells_count)=trace(corr(all_maps_first_map(cells_with_place_field(other_cells_this_bin),:)',fliplr(all_maps_second_map(cells_with_place_field(shuffle_cells_this_bin(1:length(other_cells_this_bin))),:))'))/length(shuffle_cells_this_bin);
                all_cells_stability{comparison_count}(cells_count)=corr(all_maps_first_map(cells_with_place_field(m),:)',fliplr(all_maps_second_map(cells_with_place_field(m),:))');
            end
        end
        if cells_count<length(cells_with_place_field)
            all_cells_predicted_stability{comparison_count}(cells_count+1:end)=[];
            all_cells_predicted_stability_shuffle{comparison_count}(cells_count+1:end)=[];
            all_cells_stability{comparison_count}(cells_count+1:end)=[];
        end
        
        all_pv_correlations_temp=corr(all_maps_second_map,all_maps_second_map);
        all_pv_correlations_both_sides_temp=corr(squeeze(averaged_smoothed_rate_maps_both_sides(k,:,:)),squeeze(averaged_smoothed_rate_maps_both_sides(k,:,:)));
        PV_correlations_per_map_bin_resolution((k-1)*2*length(bins_to_use)+1:k*2*length(bins_to_use),(k-1)*2*length(bins_to_use)+1:k*2*length(bins_to_use))=all_pv_correlations_temp;
        PV_correlations_per_map_both_sides((k-1)*length(bins_to_use)+1:k*length(bins_to_use),(k-1)*length(bins_to_use)+1:k*length(bins_to_use))=all_pv_correlations_both_sides_temp;
        all_pv_correlations_temp=corr(all_maps_first_map,all_maps_second_map);
        all_pv_correlations_both_sides_temp=corr(squeeze(averaged_smoothed_rate_maps_both_sides(n,:,:)),squeeze(averaged_smoothed_rate_maps_both_sides(k,:,:)));
        PV_correlations_per_map_bin_resolution((n-1)*2*length(bins_to_use)+1:n*2*length(bins_to_use),(k-1)*2*length(bins_to_use)+1:k*2*length(bins_to_use))=all_pv_correlations_temp;
        PV_correlations_per_map_bin_resolution((k-1)*2*length(bins_to_use)+1:k*2*length(bins_to_use),(n-1)*2*length(bins_to_use)+1:n*2*length(bins_to_use))=all_pv_correlations_temp';
        PV_correlations_per_map_both_sides((n-1)*length(bins_to_use)+1:n*length(bins_to_use),(k-1)*length(bins_to_use)+1:k*length(bins_to_use))=all_pv_correlations_both_sides_temp;
        PV_correlations_per_map_both_sides((k-1)*length(bins_to_use)+1:k*length(bins_to_use),(n-1)*length(bins_to_use)+1:n*length(bins_to_use))=all_pv_correlations_both_sides_temp';
        
        % calculating the rate map correlation for same versus flipped
        % orientations:
        all_rate_map_correlations_temp=corr(all_maps_first_map',all_maps_second_map');
        all_rate_map_correlations_flipped_temp=corr(fliplr(all_maps_first_map)',all_maps_second_map');
        all_concatenated_rate_map_correlations_flipped_temp=corr(fliplr(all_concatenated_maps_first_map)',all_concatenated_maps_second_map');
        
        for m=1:size(all_rate_map_correlations_temp,1)
            [this_max_value,this_max_bin]=max(all_maps_first_map(m,:));
            if this_max_value>0
                rate_map_correlations_across_maps{comparison_count,this_max_bin}=[rate_map_correlations_across_maps{comparison_count,this_max_bin},all_rate_map_correlations_temp(m,m)];
                rate_map_correlations_across_maps_flipped{comparison_count,this_max_bin}=[rate_map_correlations_across_maps_flipped{comparison_count,this_max_bin},all_rate_map_correlations_flipped_temp(m,m)];
                all_concatenated_rate_map_correlations{comparison_count}=[all_concatenated_rate_map_correlations{comparison_count},all_concatenated_rate_map_correlations_flipped_temp(m,m)];
            end
        end
    end
end

comparison_number=1;
smoothing_kernel=gausswin(round(length(bins_to_use)/4));
same_direction_right_PV_correlation=zeros(sum(1:number_of_maps-1),length(bins_to_use));
same_direction_left_PV_correlation=zeros(sum(1:number_of_maps-1),length(bins_to_use));
opposite_direction_right_left_PV_correlation=zeros(sum(1:number_of_maps-1),length(bins_to_use));
opposite_direction_left_right_PV_correlation=zeros(sum(1:number_of_maps-1),length(bins_to_use));
right_PV_correlations=zeros(sum(1:number_of_maps-1),length(bins_to_use));
left_PV_correlations=zeros(sum(1:number_of_maps-1),length(bins_to_use));
same_PV_correlation=zeros(sum(1:number_of_maps-1),length(bins_to_use));
opposite_PV_correlation=zeros(sum(1:number_of_maps-1),length(bins_to_use));
comparison_count=0;
all_PV_correlations=[];

for l=1:number_of_maps
    for m=l+1:number_of_maps
        comparison_count=comparison_count+1;
        
        % fraction of stable cells:
        correlation_threshold=0.5;
        all_start_same=[];
        all_start_flipped=[];
        for n=1:round(length(bins_to_use)/3)
            all_start_same=[all_start_same, rate_map_correlations_across_maps{comparison_count,n}];
            all_start_same=[all_start_same, rate_map_correlations_across_maps{comparison_count,2*length(bins_to_use)+1-n}];
            all_start_flipped=[all_start_flipped, rate_map_correlations_across_maps_flipped{comparison_count,n}];
            all_start_flipped=[all_start_flipped, rate_map_correlations_across_maps_flipped{comparison_count,2*length(bins_to_use)+1-n}];
        end
        fraction_stable_same_start(comparison_count)=sum(all_start_same>correlation_threshold)/length(all_start_same);
        fraction_stable_flipped_start(comparison_count)=sum(all_start_flipped>correlation_threshold)/length(all_start_flipped);
        
        all_end_same=[];
        all_end_flipped=[];
        for n=length(bins_to_use)-round(length(bins_to_use)/3)+1:length(bins_to_use)
            all_end_same=[all_end_same, rate_map_correlations_across_maps{comparison_count,n}];
            all_end_same=[all_end_same, rate_map_correlations_across_maps{comparison_count,2*length(bins_to_use)+1-n}];
            all_end_flipped=[all_end_flipped, rate_map_correlations_across_maps_flipped{comparison_count,n}];
            all_end_flipped=[all_end_flipped, rate_map_correlations_across_maps_flipped{comparison_count,2*length(bins_to_use)+1-n}];
        end
        fraction_stable_same_end=sum(all_end_same>correlation_threshold)/length(all_end_same);
        fraction_stable_flipped_end=sum(all_end_flipped>correlation_threshold)/length(all_end_flipped);
        
        all_middle_same=[];
        all_middle_flipped=[];
        for n=round(length(bins_to_use)/3)+1:length(bins_to_use)-round(length(bins_to_use)/3)
            all_middle_same=[all_middle_same, rate_map_correlations_across_maps{comparison_count,n}];
            all_middle_same=[all_middle_same, rate_map_correlations_across_maps{comparison_count,2*length(bins_to_use)+1-n}];
            all_middle_flipped=[all_middle_flipped, rate_map_correlations_across_maps_flipped{comparison_count,n}];
            all_middle_flipped=[all_middle_flipped, rate_map_correlations_across_maps_flipped{comparison_count,2*length(bins_to_use)+1-n}];
        end
        fraction_stable_same_middle=sum(all_middle_same>correlation_threshold)/length(all_middle_same);
        fraction_stable_flipped_middle=sum(all_middle_flipped>correlation_threshold)/length(all_middle_flipped);
        
        % main diagonal and reverse diagonal:
        for n=1:length(bins_to_use)
            same_direction_right_PV_correlation(comparison_count,n)=PV_correlations_per_map_bin_resolution((l-1)*2*length(bins_to_use)+n,(m-1)*2*length(bins_to_use)+n);
            same_direction_left_PV_correlation(comparison_count,n)=PV_correlations_per_map_bin_resolution((l-1)*2*length(bins_to_use)+length(bins_to_use)+n,(m-1)*2*length(bins_to_use)+length(bins_to_use)+n);
            opposite_direction_right_left_PV_correlation(comparison_count,n)=PV_correlations_per_map_bin_resolution((l-1)*2*length(bins_to_use)+length(bins_to_use)*2+1-n,(m-1)*2*length(bins_to_use)+n);
            opposite_direction_left_right_PV_correlation(comparison_count,n)=PV_correlations_per_map_bin_resolution((l-1)*2*length(bins_to_use)+length(bins_to_use)*1+1-n,(m-1)*2*length(bins_to_use)+length(bins_to_use)+n);
        end
        
        same_average_PV_correlation_first_half(comparison_count)=mean([same_direction_right_PV_correlation(comparison_count,1:round(length(same_direction_right_PV_correlation)/2)), same_direction_left_PV_correlation(comparison_count,1:round(length(same_direction_left_PV_correlation)/2))]);
        same_average_PV_correlation_second_half(comparison_count)=mean([same_direction_right_PV_correlation(comparison_count,round(length(same_direction_right_PV_correlation)/2)+1:end), same_direction_left_PV_correlation(comparison_count,round(length(same_direction_right_PV_correlation)/2)+1:end)]);
        opposite_average_PV_correlation_first_half(comparison_count)=mean([opposite_direction_right_left_PV_correlation(comparison_count,1:round(length(opposite_direction_right_left_PV_correlation)/2)), opposite_direction_left_right_PV_correlation(comparison_count,round(length(opposite_direction_left_right_PV_correlation)/2)+1:end)]);
        opposite_average_PV_correlation_second_half(comparison_count)=mean([opposite_direction_right_left_PV_correlation(comparison_count,round(length(opposite_direction_right_left_PV_correlation)/2)+1:end), opposite_direction_left_right_PV_correlation(comparison_count,1:round(length(opposite_direction_left_right_PV_correlation)/2))]);
        same_PV_correlation(comparison_count,:)=mean([same_direction_right_PV_correlation(comparison_count,:); same_direction_left_PV_correlation(comparison_count,:)]);
        opposite_PV_correlation(comparison_count,:)=mean([opposite_direction_right_left_PV_correlation(comparison_count,:); fliplr(opposite_direction_left_right_PV_correlation(comparison_count,:))]);
        
        same_std_PV_correlation_first_half(comparison_count)=std([same_direction_right_PV_correlation(comparison_count,1:round(length(same_direction_right_PV_correlation)/2)), same_direction_left_PV_correlation(comparison_count,1:round(length(same_direction_left_PV_correlation)/2))]);
        same_std_PV_correlation_second_half(comparison_count)=std([same_direction_right_PV_correlation(comparison_count,round(length(same_direction_right_PV_correlation)/2)+1:end), same_direction_left_PV_correlation(comparison_count,round(length(same_direction_right_PV_correlation)/2)+1:end)]);
        opposite_std_PV_correlation_first_half(comparison_count)=std([opposite_direction_right_left_PV_correlation(comparison_count,1:round(length(opposite_direction_right_left_PV_correlation)/2)), opposite_direction_left_right_PV_correlation(comparison_count,round(length(opposite_direction_left_right_PV_correlation)/2)+1:end)]);
        opposite_std_PV_correlation_second_half(comparison_count)=std([opposite_direction_right_left_PV_correlation(comparison_count,round(length(opposite_direction_right_left_PV_correlation)/2)+1:end), opposite_direction_left_right_PV_correlation(comparison_count,1:round(length(opposite_direction_left_right_PV_correlation)/2))]);
        same_direction_right_PV_correlation(comparison_count,:)=conv(same_direction_right_PV_correlation(comparison_count,:),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        same_direction_left_PV_correlation(comparison_count,:)=conv(same_direction_left_PV_correlation(comparison_count,:),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        opposite_direction_right_left_PV_correlation(comparison_count,:)=conv(opposite_direction_right_left_PV_correlation(comparison_count,:),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        opposite_direction_left_right_PV_correlation(comparison_count,:)=conv(opposite_direction_left_right_PV_correlation(comparison_count,:),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        
        % same position on opposite directions:
        for n=1:length(bins_to_use)
            right_PV_correlations(comparison_count,n)=PV_correlations_per_map_bin_resolution((l-1)*2*length(bins_to_use)+n,(m)*2*length(bins_to_use)+1-n);
            left_PV_correlations(comparison_count,n)=PV_correlations_per_map_bin_resolution((l-1)*2*length(bins_to_use)+length(bins_to_use)+n,(m)*2*length(bins_to_use)+1-length(bins_to_use)-n);
        end
    end
end

% comparison of rate map correlations between maps:
correlation_values=-0.35:0.1:0.95;
rate_map_correlations_in_versus_another_map=cell(1,length(correlation_values));
rate_map_relations_across_maps=[];
comparison_count=0;
for n=1:length(all_significant_both_directions)
    for k=n+1:length(all_significant_both_directions)
        for m=1:number_of_maps
            first_map_correlation=rate_map_correlations_per_map{m}(all_significant_both_directions(n),all_significant_both_directions(k));
            [~,first_map_correlation_ind]=min(abs(first_map_correlation-correlation_values));
            for l=m+1:number_of_maps
                second_map_correlation=rate_map_correlations_per_map{l}(all_significant_both_directions(n),all_significant_both_directions(k));
                if ~isnan(first_map_correlation) & ~isnan(second_map_correlation)
                    comparison_count=comparison_count+1;
                    rate_map_correlations_in_versus_another_map{first_map_correlation_ind}=[rate_map_correlations_in_versus_another_map{first_map_correlation_ind},second_map_correlation];
                    rate_map_relations_across_maps(1,comparison_count)=first_map_correlation;
                    rate_map_relations_across_maps(2,comparison_count)=second_map_correlation;
                end
            end
        end
    end
end
mean_rate_map_correlation=zeros(1,length(correlation_values));
std_rate_map_correlation=zeros(1,length(correlation_values));
number_of_values=zeros(1,length(correlation_values));
for n=1:length(correlation_values)
    mean_rate_map_correlation(n)=mean(rate_map_correlations_in_versus_another_map{n});
    std_rate_map_correlation(n)=std(rate_map_correlations_in_versus_another_map{n});
    number_of_values(n)=length(rate_map_correlations_in_versus_another_map{n});
end
figure
errorbar(correlation_values,mean_rate_map_correlation,std_rate_map_correlation./sqrt(number_of_values),'linewidth',2);
xlabel('Correlation in a given map')
ylabel('Correlation in another map')
set(gca,'fontsize',16)
box off
axis square

% plotting the bin resolution PV correlations:
figure('units','normalized','outerposition',[0.35 0.2 0.3 0.6])
axes('position',[0.16 0.2 0.7 0.7])
imagesc(PV_correlations_per_map_bin_resolution)
set(gca,'xtick',[])
set(gca,'ytick',[])
min_correlation=min(min(PV_correlations_per_map_bin_resolution));
caxis([round(min_correlation*100)/100 1])
axis square
colormap('jet')
for n=1:number_of_maps
    hold on
    plot(0.5+[0 number_of_maps*2*length(bins_to_use)],0.5+[(n*2-1)*length(bins_to_use) (n*2-1)*length(bins_to_use)],'linewidth',2,'color','k')
    hold on
    plot(0.5+[(n*2-1)*length(bins_to_use) (n*2-1)*length(bins_to_use)],0.5+[0 number_of_maps*2*length(bins_to_use)],'linewidth',2,'color','k')
end
for n=1:number_of_maps
    hold on
    if n<number_of_maps
        plot(0.5+[0 number_of_maps*2*length(bins_to_use)],0.5+[n*2*length(bins_to_use) n*2*length(bins_to_use)],'linewidth',3,'color','w')
        hold on
        plot(0.5+[n*2*length(bins_to_use) n*2*length(bins_to_use)],0.5+[0 number_of_maps*2*length(bins_to_use)],'linewidth',3,'color','w')
        hold on
    end
    text(0.5+2*n*length(bins_to_use)-length(bins_to_use),number_of_maps*2*length(bins_to_use)+15,['Map ' num2str(n)],'fontsize',16,'horizontalalignment','center')
    text(-15,0.5+2*n*length(bins_to_use)-length(bins_to_use),['Map ' num2str(n)],'fontsize',16,'rotation',90,'verticalalignment','middle','horizontalalignment','center')
    text(0.5+2*n*length(bins_to_use)-1.5*length(bins_to_use),number_of_maps*2*length(bins_to_use)+5,'Right','fontsize',16,'horizontalalignment','center')
    text(0.5+2*n*length(bins_to_use)-0.5*length(bins_to_use),number_of_maps*2*length(bins_to_use)+5,'Left','fontsize',16,'horizontalalignment','center')
    text(-5,0.5+2*n*length(bins_to_use)-1.5*length(bins_to_use),'Right','fontsize',16,'rotation',90,'verticalalignment','middle','horizontalalignment','center')
    text(-5,0.5+2*n*length(bins_to_use)-0.5*length(bins_to_use),'Left','fontsize',16,'rotation',90,'verticalalignment','middle','horizontalalignment','center')
end

set(gca,'fontsize',16)
axes('position',[0.88 0.23 0.03 0.64])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'PV correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,num2str(round(min_correlation*100)/100),'fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)
savefig(fullfile(figures_directory,'PV correlations across maps.fig'))
saveas(gcf,fullfile(figures_directory,'PV correlations across maps'),'png')

% plotting the main and reverse diagonals:
if number_of_maps>1
    x_vec=bin_size/2:bin_size:length(bins_to_use)*bin_size;
    figure
    plot(x_vec,same_direction_right_PV_correlation(comparison_number,:),'color','b','linewidth',2)
    hold on
    plot(x_vec,same_direction_left_PV_correlation(comparison_number,:),'color','c','linewidth',2)
    hold on
    plot(x_vec,opposite_direction_right_left_PV_correlation(comparison_number,:),'color','r','linewidth',2)
    hold on
    plot(x_vec,opposite_direction_left_right_PV_correlation(comparison_number,:),'color','m','linewidth',2)
    xlabel('Position (cm)')
    ylabel('PV correlation')
    legend('Same direction (right)','Same direction (left)','Opposite direction (right vs. left)','Opposite direction (left vs. right)')
    legend('boxoff')
    box off
    set(gca,'fontsize',16)
    savefig(fullfile(figures_directory,'PV correlations main diagonals.fig'))
    saveas(gcf,fullfile(figures_directory,'PV correlations main diagonals'),'png')
    
    % plotting the correlations in the begining versus ends of main/reverse:
    figure('units','normalized','outerposition',[0.35 0.3 0.3 0.4])
    axes('position',[0.1 0.2 0.8 0.75])
    bar([1 2 3 4],[same_average_PV_correlation_first_half(comparison_number), same_average_PV_correlation_second_half(comparison_number),opposite_average_PV_correlation_first_half(comparison_number),opposite_average_PV_correlation_second_half(comparison_number)],0.7,'FaceColor','none')
    hold on
    errorbar([1 2 3 4],[same_average_PV_correlation_first_half(comparison_number), same_average_PV_correlation_second_half(comparison_number),opposite_average_PV_correlation_first_half(comparison_number),opposite_average_PV_correlation_second_half(comparison_number)],[same_std_PV_correlation_first_half(comparison_number), same_std_PV_correlation_second_half(comparison_number),opposite_std_PV_correlation_first_half(comparison_number),opposite_std_PV_correlation_second_half(comparison_number)],'.','linewidth',3,'color','k')
    xlim([0 5])
    set(gca,'xticklabels',{'SS','SE','FS','FE'})
    ylabel('PV correlation')
    set(gca,'fontsize',16)
    box off
end

%% Looking for fast switches between spatial representations:
event_smoothing_kernel=ones(1,3);
minimal_number_of_active_neurons=3;
frame_correlations_per_trial=cell(number_of_maps,num_sessions,num_trials);
frame_correlations_concatenated=zeros(number_of_maps,1000000);
all_frame_correlations_other_maps=[];
all_frame_correlations_this_map=[];
frame_count=0;
for n=1:num_sessions
    for k=1:num_trials
        this_trial_map=map_number_per_trial(n,k);
        this_trial_events=binary_events_running_trials{n,k};
        this_trial_smoothed_events=zeros(size(this_trial_events));
        for l=1:number_of_neurons
            this_trial_smoothed_events(:,l)=conv(double(this_trial_events(:,l)),event_smoothing_kernel,'same');
        end
        this_trial_position=position_running_trials{n,k};
        this_trial_velocity=velocity_running_trials{n,k};
        for m=1:size(this_trial_events,1)
            this_frame_activity=this_trial_smoothed_events(m,intersect(active_cells_trials{n,k},significant_trials{n,k}));
            this_frame_position=this_trial_position(m)-min(bins_to_use)+1;
            this_frame_velocity=this_trial_velocity(m);
            this_frame_direction=-sign(this_frame_velocity)/2+1.5;
            if ~isempty(intersect(this_frame_position,bins_to_use-min(bins_to_use)+1)) & sum(this_frame_activity)>=minimal_number_of_active_neurons
                frame_count=frame_count+1;
                for p=1:number_of_maps
                    temp_frame_corr=corrcoef(this_frame_activity,squeeze(averaged_smoothed_rate_maps(p,intersect(active_cells_trials{n,k},significant_trials{n,k}),this_frame_position,this_frame_direction)));
                    frame_correlations_per_trial{p,n,k}=[frame_correlations_per_trial{p,n,k},temp_frame_corr(1,2)];
                    frame_correlations_concatenated(p,frame_count)=temp_frame_corr(1,2);
                    if this_trial_map==p
                        all_frame_correlations_this_map=[all_frame_correlations_this_map,temp_frame_corr(1,2)];
                    else
                        all_frame_correlations_other_maps=[all_frame_correlations_other_maps,temp_frame_corr(1,2)];
                    end
                end
            end
        end
    end
end
frame_correlations_concatenated(:,frame_count+1:end)=[];

% calculating within trial correlations:
chosen_session=1;
smoothing_kernel=gausswin(1+round(length(bins_to_use)/4));
time_window=60;
correlation_resolution=1; % in seconds
auto_correlation_vector=zeros(3,floor((size(concatanated_binary_events{chosen_session},1)-fs*time_window)/(fs*correlation_resolution)));
t_auto_correlation=time_window/2:correlation_resolution:time_window/2+correlation_resolution*(length(auto_correlation_vector)-1);
for n=1:length(auto_correlation_vector)
    window_indexes=(t_auto_correlation(n)-time_window/2)*fs+1:(t_auto_correlation(n)+time_window/2)*fs;
    this_comparison_number_of_neurons=length(intersect(active_cells_sessions{chosen_session},significant_sessions{chosen_session}));
    this_comparison_neurons=intersect(active_cells_sessions{chosen_session},significant_sessions{chosen_session});
    rate_maps_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
    smoothed_rate_maps_window=zeros(size(rate_maps_window));
    
    for k=1:length(bins_to_use)
        this_bins_right_index=find(round(concatanated_locations{chosen_session}(window_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(window_indexes)>velocity_thresh);
        this_bins_left_index=find(round(concatanated_locations{chosen_session}(window_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(window_indexes)<-velocity_thresh);
        rate_maps_window(:,k,1)=mean(concatanated_binary_events{chosen_session}(window_indexes(1)-1+this_bins_right_index,this_comparison_neurons));
        rate_maps_window(:,k,2)=mean(concatanated_binary_events{chosen_session}(window_indexes(1)-1+this_bins_left_index,this_comparison_neurons));
    end
    
    for m=1:this_comparison_number_of_neurons
        smoothed_rate_maps_window(m,:,1)=conv(rate_maps_window(m,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        smoothed_rate_maps_window(m,:,2)=conv(rate_maps_window(m,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
    end
    
    smoothed_rate_maps_window_both_sides=[squeeze(smoothed_rate_maps_window(:,:,1)) ; squeeze(smoothed_rate_maps_window(:,:,2))];
    
    for p=1:number_of_maps
        this_averaged_smoothed_rate_maps=[squeeze(averaged_smoothed_rate_maps(p,this_comparison_neurons,:,1)) ; squeeze(averaged_smoothed_rate_maps(p,this_comparison_neurons,:,2))];
        all_PV_correlations_temp=corr(smoothed_rate_maps_window_both_sides,this_averaged_smoothed_rate_maps);
        PV_correlations_temp=zeros(1,length(bins_to_use));
        for k=1:length(bins_to_use)
            PV_correlations_temp(k)=all_PV_correlations_temp(k,k);
        end
        auto_correlation_vector(p,n)=mean(PV_correlations_temp,'omitnan');
    end
end

% plotting the distribution of PV correlations for the same map versus
% different maps:
auto_correlation_smoothing_kernel=gausswin(9);
smoothed_auto_correlation_vector=zeros(size(t_auto_correlation));
for n=1:number_of_maps
    smoothed_auto_correlation_vector(n,:)=conv(auto_correlation_vector(n,:),auto_correlation_smoothing_kernel,'same')./conv(ones(1,length(auto_correlation_vector(n,:))),auto_correlation_smoothing_kernel,'same');
end
figure
for n=1:number_of_maps
    plot(t_auto_correlation,smoothed_auto_correlation_vector(n,:),'color',color_vector(n,:),'linewidth',2)
    hold on
end
for n=1:num_trials-1
    plot([180*n 180*n],[-0.1 1],'--','color','k','linewidth',2)
end
xlim([0 t_auto_correlation(end-7)])
ylim([-0.1 1])
xlabel('Time (sec)')
ylabel('PV correlation')
set(gca,'fontsize',16)
legend('Map 1','Map 2')
legend('boxoff')
box off

figure
for n=1:number_of_maps
    plot(frame_correlations_concatenated(n,:),'color',color_vector(n,:),'linewidth',2)
    hold on
end

figure
x_vec=-0.195:0.01:0.995;
[n1,~]=hist(all_frame_correlations_this_map,x_vec);
[n2,~]=hist(all_frame_correlations_other_maps,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
bar(x_vec,n1,'b','barwidth',1,'FaceAlpha',0.6);
hold on
bar(x_vec,n2,'r','barwidth',1,'FaceAlpha',0.6);
xlim([-0.2 1])
xlabel('PV correlation')
ylabel('Frequency (frames)')
set(gca,'fontsize',16)
legend('Same map','Different maps')
legend('boxoff')
box off

%% Map switch detector based on PV correlations to a reference map:
same_map_threshold=0.2;
different_maps_threshold=0.2;
end_trial_window_size=60;
smoothing_kernel=gausswin(1+round(num_bins/4));
time_window=60;
correlation_resolution=3; % in seconds

% calculating within trial correlations:
num_switches_within_trials=0;
num_switches_across_trials=0;
for s=1:num_sessions
    chosen_session=s;    
    auto_correlation_vector=zeros(3,floor((size(concatanated_binary_events{chosen_session},1)-fs*time_window)/(fs*correlation_resolution)));
    map_vector=zeros(1,floor((size(concatanated_binary_events{chosen_session},1)-fs*time_window)/(fs*correlation_resolution)));
    map_vector_2=zeros(1,floor((size(concatanated_binary_events{chosen_session},1)-fs*time_window)/(fs*correlation_resolution)));
    t_auto_correlation=time_window/2:correlation_resolution:time_window/2+correlation_resolution*(length(auto_correlation_vector)-1);
    
    time_of_trials=frame_log{s}(relevant_trials,3:4)/fs;
    time_of_trials_plus_window=time_of_trials+repmat([-end_trial_window_size/2 end_trial_window_size/2],num_trials,1);
    end_trial_window_range=zeros(length(relevant_trials)-1,2);
    for l=1:num_trials-1
        end_trial_window_range(l,1)=time_of_trials_plus_window(l+1,1);
        end_trial_window_range(l,2)=time_of_trials_plus_window(l,2);
    end
    end_trial_window_range=end_trial_window_range-time_of_trials(1,1);    
    map_found=0;
    for n=1:length(auto_correlation_vector)
        window_indexes=(t_auto_correlation(n)-time_window/2)*fs+1:(t_auto_correlation(n)+time_window/2)*fs;
        this_comparison_number_of_neurons=length(intersect(active_cells_sessions{chosen_session},significant_sessions{chosen_session}));
        this_comparison_neurons=intersect(active_cells_sessions{chosen_session},significant_sessions{chosen_session});
        rate_maps_window=zeros(this_comparison_number_of_neurons,length(bins_to_use),2);
        smoothed_rate_maps_window=zeros(size(rate_maps_window));
        
        for k=1:length(bins_to_use)
            this_bins_right_index=find(round(concatanated_locations{chosen_session}(window_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(window_indexes)>velocity_thresh);
            this_bins_left_index=find(round(concatanated_locations{chosen_session}(window_indexes))==bins_to_use(k) & concatanated_velocities{chosen_session}(window_indexes)<-velocity_thresh);
            rate_maps_window(:,k,1)=mean(concatanated_binary_events{chosen_session}(window_indexes(1)-1+this_bins_right_index,this_comparison_neurons),'omitnan');
            rate_maps_window(:,k,2)=mean(concatanated_binary_events{chosen_session}(window_indexes(1)-1+this_bins_left_index,this_comparison_neurons),'omitnan');
        end
        
        for m=1:this_comparison_number_of_neurons
            smoothed_rate_maps_window(m,:,1)=conv(rate_maps_window(m,:,1),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
            smoothed_rate_maps_window(m,:,2)=conv(rate_maps_window(m,:,2),smoothing_kernel,'same')./conv(ones(1,length(bins_to_use)),smoothing_kernel,'same');
        end
        
        smoothed_rate_maps_window_both_sides=[squeeze(smoothed_rate_maps_window(:,:,1)) ; squeeze(smoothed_rate_maps_window(:,:,2))];
        
        for p=1:number_of_maps
            this_averaged_smoothed_rate_maps=[squeeze(averaged_smoothed_rate_maps(p,this_comparison_neurons,:,1)) ; squeeze(averaged_smoothed_rate_maps(p,this_comparison_neurons,:,2))];
            all_PV_correlations_temp=corr(smoothed_rate_maps_window_both_sides,this_averaged_smoothed_rate_maps);
            PV_correlations_temp=zeros(1,length(bins_to_use));
            for k=1:length(bins_to_use)
                PV_correlations_temp(k)=all_PV_correlations_temp(k,k);
            end
            auto_correlation_vector(p,n)=mean(PV_correlations_temp,'omitnan');            
        end
    end
    for n=1:length(auto_correlation_vector)
        [sorted_correlations,max_map]=sort(auto_correlation_vector(:,n),'descend');
        max_correlation=sorted_correlations(1);
        second_max_correlation=sorted_correlations(2);
        if max_correlation>same_map_threshold    
            if map_found==0
                map_found=1;
            end
            map_vector(n)=max_map(1);
            map_vector_2(n)=max_map(1);
        else
            if map_found
                map_vector(n)=map_vector(n-1);
            end
        end
    end
    map_switch_vector=diff(map_vector);
    map_switch_vector(map_vector(1:end-1)==0)=0;   
    map_switch_vector(map_switch_vector~=0 & map_vector(2:end)~=0 & map_vector(1:end-1)~=0)=0.8;
    
    non_zero_map_vector=map_vector_2(map_vector_2>0);
    non_zero_t_auto_correlation=t_auto_correlation(map_vector_2>0);
    non_zero_auto_correlation_vector=auto_correlation_vector(:,map_vector_2>0);
     
    map_switch_vector_2=zeros(size(non_zero_t_auto_correlation));
    for l=1:num_trials
        this_trial_indexes=find(non_zero_t_auto_correlation>time_of_trials(l,1)+end_trial_window_size/2-time_of_trials(1,1) & non_zero_t_auto_correlation<time_of_trials(l,2)-end_trial_window_size/2-time_of_trials(1,1));
        if ~isempty(this_trial_indexes)
            this_trial_map_vector=non_zero_map_vector(this_trial_indexes);
            num_switches_within_trials=num_switches_within_trials+sum(diff(this_trial_map_vector)~=0);
            map_switch_vector_2(this_trial_indexes(diff(this_trial_map_vector)~=0))=0.8;
        end
        if l<num_trials
            this_across_trials_indexes=find(non_zero_t_auto_correlation>time_of_trials(l,2)-end_trial_window_size/2-time_of_trials(1,1) & non_zero_t_auto_correlation<time_of_trials(l+1,1)+end_trial_window_size/2-time_of_trials(1,1));
            if ~isempty(this_across_trials_indexes)
                this_across_trials_map_vector=non_zero_map_vector(this_across_trials_indexes);
                if sum(diff(this_across_trials_map_vector)~=0)>0
                    num_switches_across_trials=num_switches_across_trials+1;
                    map_switch_vector_2(this_across_trials_indexes(diff(this_across_trials_map_vector)~=0))=0.8;
                end
            end
        end
    end
end

%% Quantifying spatial information per map versus shuffle:
smoothing_kernel=gausswin(1+round(length(bins_to_use)/4));
x_vec=1:length(bins_to_use);
mean_spatial_info_per_map=cell(1,number_of_maps);
num_significant_cells_per_map=cell(1,number_of_maps);
num_active_cells_per_map=cell(1,number_of_maps);
fraction_active_cells_per_map_1=cell(1,number_of_maps);
fraction_active_cells_per_map_2=cell(1,number_of_maps);
average_number_events_per_map=cell(1,number_of_maps);
fraction_significant_cells_per_map=cell(1,number_of_maps);
fraction_significant_cells_days_without_maps=[];
mean_spatial_info_days_without_maps=[];
shuffle_mean_spatial_info_per_map=cell(1,number_of_maps);
shuffle_num_significant_cells_per_map=cell(1,number_of_maps);
shuffle_fraction_significant_cells_per_map=cell(1,number_of_maps);
preferred_locations_per_map=cell(1,number_of_maps);
for n=1:num_sessions
    this_day_maps=unique(map_number_per_trial(n,:));
    for k=1:num_trials
        this_map=map_number_per_trial(n,k);
        if this_map>0
            mean_spatial_info_per_map{this_map}=[mean_spatial_info_per_map{this_map},average_spatial_information_trials(n,k)];
            num_significant_cells_per_map{this_map}=[num_significant_cells_per_map{this_map},number_of_significant_cells_trials(n,k)];
            fraction_significant_cells_per_map{this_map}=[fraction_significant_cells_per_map{this_map},fraction_significant_trials(n,k)];
            shuffle_mean_spatial_info_per_map{this_map}=[shuffle_mean_spatial_info_per_map{this_map},shuffle_spatial_information_trials(n,k)];
            shuffle_num_significant_cells_per_map{this_map}=[shuffle_num_significant_cells_per_map{this_map},shuffle_number_of_significant_cells_trials(n,k)];
            shuffle_fraction_significant_cells_per_map{this_map}=[shuffle_fraction_significant_cells_per_map{this_map},shuffle_fraction_of_significant_cells_trials(n,k)];
            this_trial_preferred_locations=[squeeze(preferred_locations_trials(n,k,:,1)) ; squeeze(preferred_locations_trials(n,k,:,2))];
            this_distribution=hist(this_trial_preferred_locations(significant_trials{n,k}),x_vec);
            this_distribution=conv(this_distribution,smoothing_kernel,'same')./conv(ones(length(bins_to_use),1),smoothing_kernel,'same');
            preferred_locations_per_map{this_map}=[preferred_locations_per_map{this_map};this_distribution];
            num_active_cells_per_map{this_map}=[num_active_cells_per_map{this_map},length(active_cells_trials{n,k})];
            fraction_active_cells_per_map_1{this_map}=[fraction_active_cells_per_map_1{this_map},length(active_cells_trials{n,k})/neurons_per_session(n)];
            fraction_active_cells_per_map_2{this_map}=[fraction_active_cells_per_map_2{this_map},length(active_cells_trials{n,k})/number_of_neurons];
            average_number_events_per_map{this_map}=[average_number_events_per_map{this_map},mean(sum(all_events{n,k}(:,active_cells_trials{n,k})>0))];
            if length(this_day_maps)<2
                fraction_significant_cells_days_without_maps=[fraction_significant_cells_days_without_maps,fraction_significant_trials(n,k)];
                mean_spatial_info_days_without_maps=[mean_spatial_info_days_without_maps,average_spatial_information_trials(n,k)];
            end
        end
    end
end

average_dist_preferred_locations_per_map=zeros(number_of_maps,length(bins_to_use));
std_dist_preferred_locations_per_map=zeros(number_of_maps,length(bins_to_use));
for n=1:number_of_maps
    average_dist_preferred_locations_per_map(n,:)=mean(preferred_locations_per_map{n});
    sem_dist_preferred_locations_per_map(n,:)=std(preferred_locations_per_map{n})/sqrt(size(preferred_locations_per_map{n},1));
end

% plotting the distribution of spatial information per map:
figure
plot(100*fraction_significant_cells_per_map{1},mean_spatial_info_per_map{1},'.','markersize',18,'color','r')
if number_of_maps>1
    hold on
    plot(100*fraction_significant_cells_per_map{2},mean_spatial_info_per_map{2},'.','markersize',18,'color','g')
end
hold on
if number_of_maps==3
    plot(100*fraction_significant_cells_per_map{3},mean_spatial_info_per_map{3},'.','markersize',18,'color','b')
    hold on
end
plot(100*shuffle_fraction_of_significant_cells_trials(:),shuffle_spatial_information_trials(:),'.','markersize',18,'color','k')
axis square
xlim([0 70])
ylim([1 4])
set(gca,'ytick',1:4)
xlabel('% significant cells')
ylabel('Spatial information (bits/event)')
set(gca,'fontsize',16)
if number_of_maps==3
    legend('Map 1','Map 2','Map 3','Shuffle')
else
    legend('Map 1','Map 2','Shuffle')
end
legend('boxoff')
box off
savefig(fullfile(figures_directory,'Spatial information across maps.fig'))
saveas(gcf,fullfile(figures_directory,'Spatial information across maps'),'png')

% plotting the distribution of numbers of cells and events per map:
figure
plot(num_active_cells_per_map{1},average_number_events_per_map{1},'.','markersize',15,'color','r')
if number_of_maps>1
    hold on
    plot(num_active_cells_per_map{2},average_number_events_per_map{2},'.','markersize',15,'color','g')
end
hold on
if number_of_maps==3
    plot(num_active_cells_per_map{3},average_number_events_per_map{3},'.','markersize',15,'color','b')
    hold on
end
axis square
xlim([100 400])
ylim([5 10])
set(gca,'ytick',6:10)
xlabel('Number of active cells')
ylabel('Average activity (events/trial)')
set(gca,'fontsize',16)
if number_of_maps==3
    legend('Map 1','Map 2','Map 3','Shuffle')
else
    legend('Map 1','Map 2','Shuffle')
end
legend('boxoff')
box off


%% Quantifying behavior per map:
velocity_thresh=1;
fraction_stationary_time_per_map=cell(1,number_of_maps);
average_running_speed_per_map=cell(1,number_of_maps);
number_of_track_traversals_per_map=cell(1,number_of_maps);
number_of_track_traversals_trial=zeros(num_sessions,num_trials);
average_running_speed_per_trial=zeros(num_trials*num_sessions,1);

for n=1:num_sessions
    for k=1:num_trials
        this_map=map_number_per_trial(n,k);
        if this_map>0
            this_trials_velocity=all_velocities{n,k};
            this_trials_location=all_locations{n,k};
            this_trials_location(this_trials_location>num_bins/6 & this_trials_location<5*num_bins/6)=[];
            this_number_of_track_traversals=sum(abs(diff(this_trials_location))>num_bins/2);
            number_of_track_traversals_per_map{this_map}=[number_of_track_traversals_per_map{this_map},this_number_of_track_traversals];
            number_of_track_traversals_trial(n,k)=this_number_of_track_traversals;
            fraction_stationary_time_per_map{this_map}=[fraction_stationary_time_per_map{this_map},100*sum(abs(this_trials_velocity)<velocity_thresh)/length(this_trials_velocity)];
            average_running_speed_per_map{this_map}=[average_running_speed_per_map{this_map},mean(abs(this_trials_velocity(abs(this_trials_velocity)>velocity_thresh)))];
            this_index=k+(n-1)*num_trials;
            average_running_speed_per_trial(this_index)=mean(abs(this_trials_velocity(abs(this_trials_velocity)>velocity_thresh)));
        end
    end
end

if number_of_maps>1
    figure('units','normalized','outerposition',[0.4 0.35 0.28 0.4])
    axes('position',[0.1 0.1 0.2 0.8])
    if number_of_maps==2
        bar([1 2],[mean(fraction_stationary_time_per_map{1}), mean(fraction_stationary_time_per_map{2})],0.7,'FaceColor','none')
        hold on
        errorbar([1 2],[mean(fraction_stationary_time_per_map{1}), mean(fraction_stationary_time_per_map{2})],[std(fraction_stationary_time_per_map{1})/sqrt(length(fraction_stationary_time_per_map{1})), std(fraction_stationary_time_per_map{2})/sqrt(length(fraction_stationary_time_per_map{2}))],'.','linewidth',3,'color','k')
        xlim([0.4 2.6])
        set(gca,'xticklabels',{'M1','M2'})
    else
        bar([1 2 3],[mean(fraction_stationary_time_per_map{1}), mean(fraction_stationary_time_per_map{2}),mean(fraction_stationary_time_per_map{3})],0.7,'FaceColor','none')
        hold on
        errorbar([1 2 3],[mean(fraction_stationary_time_per_map{1}), mean(fraction_stationary_time_per_map{2}),mean(fraction_stationary_time_per_map{3})],[std(fraction_stationary_time_per_map{1})/sqrt(length(fraction_stationary_time_per_map{1})), std(fraction_stationary_time_per_map{2})/sqrt(length(fraction_stationary_time_per_map{2})), std(fraction_stationary_time_per_map{3})/sqrt(length(fraction_stationary_time_per_map{3}))],'.','linewidth',3,'color','k')
        xlim([0.4 3.6])
        set(gca,'xticklabels',{'M1','M2','M3'})
    end
    box off
    ylabel('% stationary periods')
    set(gca,'fontsize',14)
    
    axes('position',[0.45 0.1 0.2 0.8])
    if number_of_maps==2
        bar([1 2],[mean(average_running_speed_per_map{1}), mean(average_running_speed_per_map{2})],0.7,'FaceColor','none')
        hold on
        errorbar([1 2],[mean(average_running_speed_per_map{1}), mean(average_running_speed_per_map{2})],[std(average_running_speed_per_map{1})/sqrt(length(average_running_speed_per_map{1})), std(average_running_speed_per_map{2})/sqrt(length(average_running_speed_per_map{2}))],'.','linewidth',3,'color','k')
        xlim([0.4 2.6])
        set(gca,'xticklabels',{'M1','M2'})
    else
        bar([1 2 3],[mean(average_running_speed_per_map{1}), mean(average_running_speed_per_map{2}),mean(average_running_speed_per_map{3})],0.7,'FaceColor','none')
        hold on
        errorbar([1 2 3],[mean(average_running_speed_per_map{1}), mean(average_running_speed_per_map{2}),mean(average_running_speed_per_map{3})],[std(average_running_speed_per_map{1})/sqrt(length(average_running_speed_per_map{1})), std(average_running_speed_per_map{2})/sqrt(length(average_running_speed_per_map{2})),std(average_running_speed_per_map{3})/sqrt(length(average_running_speed_per_map{3}))],'.','linewidth',3,'color','k')
        xlim([0.4 3.6])
        set(gca,'xticklabels',{'M1','M2','M3'})
    end
    box off
    ylabel('Average running speed (cm/sec)')
    set(gca,'fontsize',14)
    
    axes('position',[0.79 0.1 0.2 0.8])
    if number_of_maps==2
        bar([1 2],[mean(number_of_track_traversals_per_map{1}), mean(number_of_track_traversals_per_map{2})],0.7,'FaceColor','none')
        hold on
        errorbar([1 2],[mean(number_of_track_traversals_per_map{1}), mean(number_of_track_traversals_per_map{2})],[std(number_of_track_traversals_per_map{1})/sqrt(length(number_of_track_traversals_per_map{1})), std(number_of_track_traversals_per_map{2})/sqrt(length(number_of_track_traversals_per_map{2}))],'.','linewidth',3,'color','k')
        xlim([0.4 2.6])
        set(gca,'xticklabels',{'M1','M2'})
    else
        bar([1 2 3],[mean(number_of_track_traversals_per_map{1}), mean(number_of_track_traversals_per_map{2}), mean(number_of_track_traversals_per_map{3})],0.7,'FaceColor','none')
        hold on
        errorbar([1 2 3],[mean(number_of_track_traversals_per_map{1}), mean(number_of_track_traversals_per_map{2}), mean(number_of_track_traversals_per_map{3})],[std(number_of_track_traversals_per_map{1})/sqrt(length(number_of_track_traversals_per_map{1})), std(number_of_track_traversals_per_map{2})/sqrt(length(number_of_track_traversals_per_map{2})),std(number_of_track_traversals_per_map{3})/sqrt(length(number_of_track_traversals_per_map{3}))],'.','linewidth',3,'color','k')
        xlim([0.4 3.6])
        set(gca,'xticklabels',{'M1','M2','M3'})
    end
    box off
    ylabel('Number of track traversals')
    set(gca,'fontsize',14)
end
%% map occurence temporal dynamics:
% long-term dynamics (across days):
if number_of_maps>1
    maps_occurence=zeros(number_of_maps,num_sessions);
    for n=1:num_sessions
        maps_occurence(1,n)=sum(map_number_per_trial(n,:)==1);
        maps_occurence(2,n)=sum(map_number_per_trial(n,:)==2);
        if number_of_maps>2
            maps_occurence(3,n)=sum(map_number_per_trial(n,:)==3);
        end
        if number_of_maps>3
            maps_occurence(4,n)=sum(map_number_per_trial(n,:)==4);
        end
    end
    
    maps_occurence_auto_correlation=zeros(number_of_maps,2*num_sessions-1);
    for n=1:number_of_maps
        maps_occurence_auto_correlation(n,:)=normxcorr2(maps_occurence(n,:),maps_occurence(n,:));
    end
    
    average_auto_correlation=zeros(1,2*num_sessions-1);
    std_auto_correlation=zeros(1,2*num_sessions-1);
    for n=1:2*num_sessions-1
        average_auto_correlation(n)=mean(maps_occurence_auto_correlation(:,n));
        std_auto_correlation(n)=std(maps_occurence_auto_correlation(:,n));
    end
    
    % short-term dynamics (across trials):
    number_of_consecutive_data=0;
    session_count=0;
    for n=1:num_sessions
        this_maps_vector=map_number_per_trial(n,:);
        if length(unique(this_maps_vector))>1
            session_count=session_count+1;
            number_of_consecutive_data=number_of_consecutive_data+sum(diff(this_maps_vector)==0);
        end
    end
    p_consecutive_data=number_of_consecutive_data/session_count/(num_trials-1);
    
    num_shuffles=1000;
    p_consecutive_shuffle=zeros(1,num_shuffles);
    for k=1:num_shuffles
        number_of_consecutive_data=0;
        session_count=0;
        for n=1:num_sessions
            this_maps_vector=map_number_per_trial(n,:);
            [~,this_shuffle_index]=sort(rand(1,num_trials));
            this_shuffled_maps=this_maps_vector(this_shuffle_index);
            if length(unique(this_maps_vector))>1
                session_count=session_count+1;
                number_of_consecutive_data=number_of_consecutive_data+sum(diff(this_shuffled_maps)==0);
            end
        end
        p_consecutive_shuffle(k)=number_of_consecutive_data/session_count/(num_trials-1);
    end
    mean_p_consecutive_shuffle=mean(p_consecutive_shuffle);       
end

disp('Finished analyzing data set')
