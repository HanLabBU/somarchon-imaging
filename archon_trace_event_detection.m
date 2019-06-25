function archon_trace_event_detection(up_threshold, down_threshold)

    if nargin<2 || isempty(down_threshold) 
        down_threshold = 2;
    end
    
    if nargin<1 || isempty(up_threshold) 
        up_threshold = 3;
    end

    event_parameter.moving_window = 51;
    event_parameter.pre_peak_data_point = 1;
    event_parameter.post_peak_data_point = 1;
    event_parameter.noise_threshold = 3;
    event_parameter.noise_pre_extension = 0; % data points
    event_parameter.noise_post_extension = 3; % data points
    event_parameter.noise_extension = 3; % data points
    
    event_parameter.down_threshold = down_threshold;
    event_parameter.up_threshold = up_threshold;

    [selected_files,selected_folder] = uigetfile('*.mat','MultiSelect','on');
    
    if class(selected_files)=='char'
        file_list(1).name = selected_files;
    else
        file_list = cell2struct(selected_files,'name',1);
    end

    whole_tic = tic;
    
    cd(selected_folder);
    
    for file_idx=1:length(file_list)
        

        filename = file_list(file_idx).name;
        load(filename);
        
        try
            result = rmfield(result,'roi');
            result = rmfield(result,'roi_phase');
        end
        
        try
            traces = result.traces;
        catch
            fprintf(['!Skipped: ',filename,'\n']);
            continue;
        end
        
        try
            trace_time = result.trace_time;
        catch
            result.trace_time = result.camera_frame_time(1:size(result.traces,1));
            trace_time = result.trace_time;
        end
        
        for roi_idx=1:size(traces,2)
            clear event;
            event.idx=[];
            trace = traces(:,roi_idx);
            
            d_trace = diff(trace);
            d_trace = [0;d_trace];
 
            noise_idx_list = find(d_trace<(mean(d_trace)-event_parameter.noise_threshold*std(d_trace)));
            
            % connect noise index
            d_noise_idx_list = diff(noise_idx_list);
            event_parameter.noise_extension_idx = find(d_noise_idx_list<event_parameter.noise_extension);
            if ~isempty(event_parameter.noise_extension_idx)
                for idx=1:numel(event_parameter.noise_extension_idx)
                    current_idx = event_parameter.noise_extension_idx(idx);
                    noise_idx_list = [noise_idx_list;[noise_idx_list(current_idx):noise_idx_list(current_idx+1)]'];
                end
            end
            
            noise_idx_list = unique(noise_idx_list);
            
            nan_trace = trace;
            for noise_idx=noise_idx_list'
                if noise_idx-event_parameter.noise_pre_extension>0 & d_trace(noise_idx-1)<(mean(d_trace)+event_parameter.noise_threshold*std(d_trace))
                    nan_trace(noise_idx-event_parameter.noise_pre_extension:min(noise_idx+event_parameter.noise_post_extension,length(trace))) = nan;
                end
            end
            
            d_nan_trace = diff(nan_trace);
            d_nan_trace = [0;d_nan_trace];
            
            event_parameter.up_threshold_value = event_parameter.up_threshold*nanstd(d_nan_trace);
            event.event_parameter.up_threshold_value = event_parameter.up_threshold_value;
            event_parameter.down_threshold_value = event_parameter.down_threshold*nanstd(d_nan_trace);
            event.event_parameter.down_threshold_value = event_parameter.down_threshold_value;
            
            event.trace = trace;
            event.nan_trace = nan_trace;
            
            pre_d_trace = d_trace;
            if event_parameter.pre_peak_data_point>0
                for idx=1:event_parameter.pre_peak_data_point
                    shifted_d_trace = [zeros(idx,1);d_trace(1:end-idx)];
                    shifted_d_trace(shifted_d_trace<0) = 0;
                    pre_d_trace = pre_d_trace+shifted_d_trace;
                end
            end
            
            post_d_trace = d_trace;
            if event_parameter.post_peak_data_point>0
                for idx=1:event_parameter.post_peak_data_point
                    shifted_d_trace = [d_trace(idx+1:end);zeros(idx,1)];
                    shifted_d_trace(shifted_d_trace>0) = 0;
                    post_d_trace = post_d_trace+shifted_d_trace;
                end
            end
            
            up_idx_list = find(pre_d_trace>nanmean(d_nan_trace)+event_parameter.up_threshold_value);
            
            for up_idx=up_idx_list'
                if (up_idx+1)<=numel(d_trace) & d_trace(up_idx)>0 & d_trace(up_idx+1)<0 & post_d_trace(up_idx+1)<(nanmean(d_trace)-event_parameter.down_threshold_value) & ~isnan(nan_trace(up_idx))
                    event.idx = cat(1,event.idx,up_idx);
                end
            end
            
            event.time = trace_time(event.idx);
            event.roaster = zeros(size(trace));
            event.roaster(event.idx) = 1;
            
            result.roi(roi_idx).event = event;
            
        end
        
        result.event_parameter = event_parameter;

        save(filename,'result');
        fprintf(['Saved: ',filename,'\n']);
        
    end
    
    fprintf(['Total time: ',num2str(toc(whole_tic)),' seconds.\n']);

end