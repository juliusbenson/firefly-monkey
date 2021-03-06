%% add lfps
function AddLfps(this,prs)
cd(prs.filepath_neur); % cd(fullfile(prs.filepath_neur,'Sorted'));
% determine type of electrode
linearprobe_type = []; utaharray_type = [];
for k=1:length(prs.electrode_type)
    linearprobe_type = [linearprobe_type find(cellfun(@(electrode_type) strcmp(prs.electrode_type{k},electrode_type), prs.linearprobe.types),1)];
    utaharray_type = [utaharray_type find(cellfun(@(electrode_type) strcmp(prs.electrode_type{k},electrode_type), prs.utaharray.types),1)];
end

if isfield(prs, 'isRipple')
    if prs.isRipple
        
        cd(prs.filepath_neur); % cd(fullfile(prs.filepath_neur,'Sorted'));
        file_nev=dir('*.nev'); file_ns1=dir('*.ns2'); prs.neur_filetype = 'nev'; % need to check the sampling here..
        fprintf(['... reading events from ' file_nev.name '\n']);
        [events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK
        t_end_nev = [events_nev.t_end];
        t_beg_nev = [events_nev.t_beg];
        if length(t_beg_nev) > length(t_end_nev)
            if all(t_end_nev - t_beg_nev(1:length(t_end_nev)) > 0)
                events_nev.t_beg = t_beg_nev(1:length(t_end_nev));
            elseif all(t_end_nev - t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end) > 0)
                events_nev.t_beg = t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end);
            end
        end
        
        if length(this.behaviours.trials)~=length(events_nev.t_end)
            events_nev = FixEvents_nev(events_nev,this.behaviours.trials);
        end
        if length(this.behaviours.trials)==length(events_nev.t_end)
            % all seems to be good, load LFP data
            NS1 = openNSx(file_ns1.name,'read', 'uV');
            
            [~ ,electrode_id, brain_area, channels_per_area, electrode_type] = MapChannel2Electrode_Ripple(NS1.MetaTags, prs);
            
            for j=1:sum(channels_per_area)
                if NS1.MetaTags.ChannelID(1) > 1
                channel_id = NS1.MetaTags.ChannelID(j) - NS1.MetaTags.ChannelID(1) +1;
                else
                channel_id = NS1.MetaTags.ChannelID(j);
                end
                fprintf(['Segmenting LFP :: channel ' num2str(channel_id) '\n']);
                this.lfps(end+1) = lfp(channel_id,electrode_id(j),electrode_type{j});
                this.lfps(end).brain_area = brain_area{j};
                tStart = tic;
                this.lfps(end).AddTrials(NS1.Data(j,:),NS1.MetaTags.SamplingFreq,events_nev,this.behaviours,prs);
                toc(tStart)
            end
        else
            fprintf('Cannot segment LFP: Trial counts in smr and nev files do not match \n');
            fprintf(['Trial end events: NEV file - ' num2str(length(events_nev.t_end)) ...
                ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
            fprintf('Debug and try again! \n');
        end
        
    else
        if ~isempty(linearprobe_type) % assume linearprobe is recorded using Plexon
            file_plx = dir([prs.filepath_neur, '*.plx']);
            if isempty(file_plx)
               try
                   file_plx = dir([prs.filepath_neur, 'PLEXON FILES\', '*.plx']);
               catch
               end
            end
            try
                cd(fullfile(prs.filepath_neur,'PLEXON FILES','Sorted'));
            catch
                cd(fullfile(prs.filepath_neur,'Sorted'));
            end
            brain_area = prs.area{strcmp(prs.electrode_type,prs.linearprobe.types{linearprobe_type})};
            file_ead=dir('*_ead.plx'); file_lfp=dir('*_lfp.plx'); prs.neur_filetype = 'plx';
            % read events
            fprintf(['... reading events from ' file_ead.name '\n']);
            %[events_plx, prs.fs_spk] = GetEvents_plx(file_ead.name);
            if ~isempty(file_ead)
                [events_plx] = GetEvents_plx(file_ead.name);
            else
                [events_plx] = GetEvents_plx(fullfile(prs.filepath_neur, file_plx.name));
            end
            % read lfp
            %% JP: trying to see if we can salvage from replay 
            for ii = 1:size(this.behaviours.trials, 2)
                beh_t_beg(ii) = this.behaviours.trials(ii).events.t_beg;
                
            end
            
            
            
            %%
            if length(this.behaviours.trials)==length(events_plx.t_end)
                fprintf(['... reading ' file_lfp.name '\n']);
                [ch_id,electrode_id] = MapChannel2Electrode('linearprobe');
                for j=1:prs.linearprobe.channelcount(linearprobe_type)
                    fprintf(['...... channel ' num2str(j) '/' num2str(prs.linearprobe.channelcount(linearprobe_type)) '\n']);
                    if ~isempty(file_lfp)
                        [adfreq, n, ~, fn, ad] = plx_ad_v(file_lfp.name, j-1);
                    else
                        %adfreq is 20000
                        %n is 115520845 in this example...
                        %fn is 115520845 in this example...
                        %ad is the signal
                        file_lfp_name = fullfile(file_plx.folder, file_plx.name);
                        %use plx_ad_info(file_lfp_name) to make sure you are
                        %grabing the correct files..
                        %[adfreq, n, ~, fn, ad] = plx_ad_v(file_lfp_name, prs.linearprobe.channelcount(linearprobe_type)+j-1);
                        [adfreq, n, ~, fn, ad] = plx_ad_v(file_lfp_name, 32+j-1);
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % ORIGINAL CODE
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         filter_bandpass = designfilt('lowpassfir', ...        % Response type
%                             'FilterOrder',400, ...            % Filter order
%                             'PassbandFrequency',100, ...     % Frequency constraints
%                             'StopbandFrequency',150, ...
%                             'DesignMethod','ls', ...         % Design method
%                             'PassbandWeight',1, ...          % Design method options
%                             'StopbandWeight',2, ...
%                             'SampleRate',2000);               % Sample rate
%                         
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % CHANGE PARAMS
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        cut_freq = 400;
                        [b,a] = butter(4, cut_freq/(adfreq/2), 'low');
                        ad_filter = filtfilt(b, a, ad');
                        ad = ad_filter';
                    end
                    if n == fn
                        if adfreq > prs.fs_lfp
                            % 4,10 = (adfreq/prs.fs_lfp); this for downsampling of 20000 to 500
                            ad = decimate(ad,10); ad = decimate(ad,4);
                        end
                        channel_id = j;
                        fprintf(['Segmenting LFP :: channel ' num2str(channel_id) '\n']);
                        this.lfps(end+1) = lfp(channel_id,electrode_id(ch_id == channel_id),prs.linearprobe.types{linearprobe_type});
                        this.lfps(end).brain_area = brain_area;
                        this.lfps(end).AddTrials(ad,prs.fs_lfp,events_plx,this.behaviours,prs);
                    else
                        fprintf('...... LFP is fragmented. Use a machine with more RAM or contact KL\n');
                    end
                end
            else
                fprintf('Cannot segment LFP: Trial counts in smr and plx files do not match \n');
                fprintf(['Trial end events: PLX file - ' num2str(length(events_plx.t_end)) ...
                    ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
                fprintf('Debug and try again! \n');
            end
        else
            fprintf('No Plexon neural data files in the specified path \n');
        end
        
        
        if ~isempty(utaharray_type) % assume utaharray is recorded using Cereplex or Ripple
            cd(fullfile(prs.filepath_neur,'Sorted'));
            if prs.isRipple
                file_nev=dir('*.nev'); file_ns1=dir('*.ns5'); prs.neur_filetype = 'nev';
            else
                file_nev=dir('*.nev'); file_ns1=dir('*.ns1'); prs.neur_filetype = 'nev';
            end
            fprintf(['... reading events from ' file_nev.name '\n']);
            [events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK
            t_end_nev = [events_nev.t_end];
            t_beg_nev = [events_nev.t_beg];
            if length(t_beg_nev) > length(t_end_nev)
                if all(t_end_nev - t_beg_nev(1:length(t_end_nev)) > 0)
                    events_nev.t_beg = t_beg_nev(1:length(t_end_nev));
                elseif all(t_end_nev - t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end) > 0)
                    events_nev.t_beg = t_beg_nev(length(t_beg_nev)-length(t_end_nev)+1:end);
                end
            end
            
            if length(this.behaviours.trials)~=length(events_nev.t_end)
                events_nev = FixEvents_nev(events_nev,this.behaviours.trials);
            end
            if length(this.behaviours.trials)==length(events_nev.t_end)
                NS1 = openNSx(file_ns1.name,'read', 'uV');
                if NS1.MetaTags.ChannelCount ~= prs.utaharray.channelcount(utaharray_type)
                    warning('Unexpected channel count in the file \n');
                end
                if prs.isRipple
                    [ch_id,electrode_id] = MapChannel2Electrode_Ripple(NS1.MetaTags, prs);
                else
                    [ch_id,electrode_id] = MapChannel2Electrode(prs.utaharray.types{utaharray_type});
                end
                brain_area = prs.area{strcmp(prs.electrode_type,prs.utaharray.types{utaharray_type})};
                
                for j=1:prs.utaharray.channelcount(utaharray_type)
                    channel_id = NS1.MetaTags.ChannelID(j);
                    fprintf(['Segmenting LFP :: channel ' num2str(channel_id) '\n']);
                    this.lfps(end+1) = lfp(channel_id,electrode_id(ch_id == channel_id),prs.utaharray.types{utaharray_type});
                    
                    if strcmp(prs.utaharray.types{utaharray_type},'utah96') 
                        this.lfps(end).brain_area = brain_area;
                    else
                        this.lfps(end).brain_area = prs.MapDualArray2BrainArea(brain_area, this.lfps(end).electrode_id); 
                    end
                    this.lfps(end).AddTrials(NS1.Data(j,:),NS1.MetaTags.SamplingFreq,events_nev,this.behaviours,prs);
                end
            else
                fprintf('Cannot segment LFP: Trial counts in smr and nev files do not match \n');
                fprintf(['Trial end events: NEV file - ' num2str(length(events_nev.t_end)) ...
                    ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
                fprintf('Debug and try again! \n');
            end
        else
            fprintf('No Cereplex neural data files in the specified path \n');
        end
        
    end
end