%input striatum matlab files (include movement, filtered movement, spiking,
%etc from Striatal videos)
filename_list1= findNestedFiles(); % insert the data file names here
filename_list2= findNestedFiles(); % insert the data file names here
filename_list3= findNestedFiles(); % insert the data file names here
filename_list4= findNestedFiles(); % insert the data file names here
filename_list5= findNestedFiles(); % insert the data file names here
filename_list6= findNestedFiles(); % insert the data file names here
filename_list7= findNestedFiles(); % insert the data file names here
filename_list8= findNestedFiles(); % insert the data file names here
filename_list9= findNestedFiles(); % insert the data file names here
filename_list10= findNestedFiles(); % insert the data file names here

cd() % switch to the folder containing data


%Calculate spikes in each 0.5 second window and average speed in each 0.5
%second window for all trials

for k=1:14
Striatum(k).motion=[];
Striatum(k).spikes=[];
end

for filename_idx = 1:numel(filename_list2)
        
         filename = filename_list2{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(1).motion=[Striatum(1).motion; motion];
         Striatum(1).spikes=[Striatum(1).spikes; spikes];
            end
end




for filename_idx = 1:numel(filename_list3)
        
         filename = filename_list3{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(2).motion=[Striatum(2).motion; motion];
         Striatum(2).spikes=[Striatum(2).spikes; spikes];
            end
end



for filename_idx = 1:numel(filename_list4)
        
         filename = filename_list4{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes1=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            spikes2=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,2));
            
            spikes3=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,3));
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(3).motion=[Striatum(3).motion; motion];
         Striatum(3).spikes=[Striatum(3).spikes; spikes1];
         Striatum(4).motion=[Striatum(4).motion; motion];
         Striatum(4).spikes=[Striatum(4).spikes; spikes2];
         Striatum(5).motion=[Striatum(5).motion; motion];
         Striatum(5).spikes=[Striatum(5).spikes; spikes3];
            end
end


for filename_idx = 1:numel(filename_list5)
        
         filename = filename_list5{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(6).motion=[Striatum(6).motion; motion];
         Striatum(6).spikes=[Striatum(6).spikes; spikes];
            end
end


for filename_idx = 1:numel(filename_list6)
        
         filename = filename_list6{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes1=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            spikes2=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,2));
                        
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(7).motion=[Striatum(7).motion; motion];
         Striatum(7).spikes=[Striatum(7).spikes; spikes1];
         Striatum(8).motion=[Striatum(8).motion; motion];
         Striatum(8).spikes=[Striatum(8).spikes; spikes2];
         
            end
end



for filename_idx = 1:numel(filename_list7)
        
         filename = filename_list7{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(9).motion=[Striatum(9).motion; motion];
         Striatum(9).spikes=[Striatum(9).spikes; spikes];
            end
end




for filename_idx = 1:numel(filename_list9)
        
         filename = filename_list9{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes1=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            spikes2=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,2));
            
            spikes3=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,3));
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(10).motion=[Striatum(10).motion; motion];
         Striatum(10).spikes=[Striatum(10).spikes; spikes1];
         Striatum(11).motion=[Striatum(11).motion; motion];
         Striatum(11).spikes=[Striatum(11).spikes; spikes2];
         Striatum(12).motion=[Striatum(12).motion; motion];
         Striatum(12).spikes=[Striatum(12).spikes; spikes3];
            end
end



for filename_idx = 1:numel(filename_list10)
        
         filename = filename_list10{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        load(filename)

        time=result.trace_time-result.trial_start_time;
        
        total_time=max(time);
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(total_time/time(window_frames))
            spikes=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(13).motion=[Striatum(13).motion; motion];
         Striatum(13).spikes=[Striatum(13).spikes; spikes];
            end
end

load('401754_aFOV1_5.mat')
time=result.camera_frame_time-result.trial_start_time;
        
        total_time=max(time);
        time1=time(10498);
        time2=time(44999)-time(18000);
        
        window=0.5; %s
        
        window_frames=round(window/(1/result.camera_frequency));
        
        
            for i=1:floor(time1/time(window_frames))
            spikes=sum(result.KP_raster(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(result.f_movement_speed_interp(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(14).motion=[Striatum(14).motion; motion];
         Striatum(14).spikes=[Striatum(14).spikes; spikes];
            end
            
            new_movement=result.f_movement_speed_interp(18000:44999);
            for i=1:floor(time2/time(window_frames))
            spikes=sum(result.KP_raster2(i*window_frames-(window_frames-1):i*window_frames,1));
            motion=mean(new_movement(i*window_frames-(window_frames-1):i*window_frames));
           
            
            %need to get spike and motion data from for loop and keep
        %accumulating it for each trial (make struct below)
            Striatum(14).motion=[Striatum(14).motion; motion];
         Striatum(14).spikes=[Striatum(14).spikes; spikes];
            end
%% calculate average firing rate per cell (supplemental fig)
Motion=[];
for i=1:length(Striatum)
Firing_Rate(i)=2*mean(Striatum(i).spikes);
Motion=[Motion; Striatum(i).motion];
end

max(Firing_Rate)
min(Firing_Rate)

%motion histograms (supplemental fig)
histogram(Motion,50, 'FaceColor', 'k')

figure;
subplot(3,3,1)
histogram(Striatum(10).motion,30, 'FaceColor', 'k')
title('FOV 1')

subplot(3,3,2)
histogram(Striatum(6).motion,30, 'FaceColor', 'k')
title('FOV 2')

subplot(3,3,3)
histogram(Striatum(2).motion,30, 'FaceColor', 'k')
title('FOV 3')

subplot(3,3,4)
histogram(Striatum(3).motion,30, 'FaceColor', 'k')
title('FOV 4')

subplot(3,3,5)
histogram(Striatum(1).motion,30,'FaceColor', 'k')
title('FOV 5')

subplot(3,3,6)
histogram(Striatum(7).motion,30, 'FaceColor', 'k')
title('FOV 6')

subplot(3,3,7)
histogram(Striatum(9).motion,30, 'FaceColor', 'k')
title('FOV 7')

subplot(3,3,8)
histogram(Striatum(13).motion,30, 'FaceColor', 'k')
title('FOV 8')

subplot(3,3,9)
histogram(Striatum(14).motion,30, 'FaceColor', 'k')
title('FOV 9')

%%
%calculate low and high motion in cm/s and statistics
high_thresh=10;
low_thresh=5;
for j=1:length(Striatum)
    Striatum(j).high_index=find(Striatum(j).motion>=high_thresh);
    Striatum(j).low_index=find(Striatum(j).motion<=low_thresh);
    Striatum(j).high_motion=Striatum(j).motion(Striatum(j).high_index);
    Striatum(j).high_m_spikes=Striatum(j).spikes(Striatum(j).high_index);
    Striatum(j).low_motion=Striatum(j).motion(Striatum(j).low_index);
    Striatum(j).low_m_spikes=Striatum(j).spikes(Striatum(j).low_index);
    Striatum(j).p_value=ranksum(Striatum(j).low_m_spikes,Striatum(j).high_m_spikes);
end

  

%% histogram of all motion (supp figure)

All_motion=[];
j=[1 2 3 6 7 9 10 13]; %only use data from each FOV, don't need each cell (if multiple cells in FOV, motion data is same)
for k=1:length(j)
    All_motion=[All_motion; Striatum(j(k)).motion];
end

figure
histogram(All_motion,50)
xlabel('Speed (cm/s)')
ylabel('Counts')

figure;
for i=1:length(j)
    subplot(3,3,i)
    histogram(Striatum(j(i)).motion,20)
end