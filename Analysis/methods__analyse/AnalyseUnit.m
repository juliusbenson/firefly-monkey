function stats = AnalyseUnit(trials_spks,trials_behv,behv_stats,lfps,prs)

%% load analysis params
x0 = prs.x0; y0 = prs.y0; % position of the subject at trial onset
dt = prs.dt; % sampling resolution (s)
temporal_binwidth = prs.temporal_binwidth;
corr_lag = prs.corr_lag;
duration_zeropad = prs.duration_zeropad;
nbootstraps = prs.nbootstraps;
peaktimewindow = prs.peaktimewindow;
minpeakprominence = prs.minpeakprominence.neural;
mintrialsforstats = prs.mintrialsforstats;
evaluate_peaks = 1; %prs.evaluate_peaks;
compute_tuning = 1; %prs.compute_tuning;
fitGAM_tuning = prs.fitGAM_tuning;
GAM_varexp = prs.GAM_varexp;
fitNNM = prs.fitNNM;
analyse_spikeLFPrelation = prs.analyse_spikeLFPrelation;
analyse_spikeLFPrelation_allLFPs = prs.analyse_spikeLFPrelation_allLFPs;
sta_window = prs.sta_window;
duration_nanpad = prs.duration_nanpad;
phase_slidingwindow = prs.phase_slidingwindow;
analyse_temporalphase = prs.analyse_temporalphase;
ntrls = length(trials_spks);

% prepare filter
filtwidth = 5; %prs.neuralfiltwidth;
t = linspace(-2*filtwidth,2*filtwidth,4*filtwidth + 1); h = exp(-t.^2/(2*filtwidth^2)); h = h/sum(h);

%% load cases
trialtypes = fields(behv_stats.trialtype);
events = cell2mat({trials_behv.events});
continuous = cell2mat({trials_behv.continuous});

%% event-aligned, trial-averaged firing rates
if evaluate_peaks
    gettuning = [prs.tuning_events {'saccade'}];
    for i=1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        disp(['event-aligned firing rates, trialtype: ' trialtypes{i}]);
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).events = stats.trialtype.all.events;
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                trials_spks_temp = trials_spks(trlindx);
                stats.trialtype.(trialtypes{i})(j).val = behv_stats.trialtype.(trialtypes{i})(j).val;
                %% aligned to movement onset
                if any(strcmp(gettuning,'move'))
                    trials_spks_temp2 = ShiftSpikes(trials_spks_temp,[events_temp.t_move]);
                    [nspk,ts] = Spiketimes2Rate(trials_spks_temp2,prs.ts.move,temporal_binwidth);
                    stats.trialtype.(trialtypes{i})(j).events.move.rate = nspk;
                    stats.trialtype.(trialtypes{i})(j).events.move.time = ts;
                    stats.trialtype.(trialtypes{i})(j).events.move.peakresp = ...           % significance of peak response
                        EvaluatePeakresponse(trials_spks_temp2,prs.ts.move,temporal_binwidth,peaktimewindow,minpeakprominence,nbootstraps,mintrialsforstats);
                end
                %% aligned to target onset
                if any(strcmp(gettuning,'target'))
                    trials_spks_temp2 = ShiftSpikes(trials_spks_temp,[events_temp.t_beg]-[events_temp.t_beg]);
                    [nspk,ts] = Spiketimes2Rate(trials_spks_temp2,prs.ts.target,temporal_binwidth);
                    stats.trialtype.(trialtypes{i})(j).events.target.rate = nspk;
                    stats.trialtype.(trialtypes{i})(j).events.target.time = ts;
                    stats.trialtype.(trialtypes{i})(j).events.target.peakresp = ...         % significance of peak response
                        EvaluatePeakresponse(trials_spks_temp2,prs.ts.target,temporal_binwidth,peaktimewindow,minpeakprominence,nbootstraps,mintrialsforstats);
                end
                %% aligned to movement stop
                if any(strcmp(gettuning,'stop'))
                    trials_spks_temp2 = ShiftSpikes(trials_spks_temp,[events_temp.t_stop]);
                    [nspk,ts] = Spiketimes2Rate(trials_spks_temp2,prs.ts.stop,temporal_binwidth);
                    stats.trialtype.(trialtypes{i})(j).events.stop.rate = nspk;
                    stats.trialtype.(trialtypes{i})(j).events.stop.time = ts;
                    stats.trialtype.(trialtypes{i})(j).events.stop.peakresp = ...           % significance of peak response
                        EvaluatePeakresponse(trials_spks_temp2,prs.ts.stop,temporal_binwidth,peaktimewindow,minpeakprominence,nbootstraps,mintrialsforstats);
                end
                %% aligned to reward
                if any(strcmp(gettuning,'reward'))
                    trials_spks_temp2 = ShiftSpikes(trials_spks_temp,[events_temp.t_rew]);
                    [nspk,ts] = Spiketimes2Rate(trials_spks_temp2,prs.ts.reward,temporal_binwidth);
                    stats.trialtype.(trialtypes{i})(j).events.reward.rate = nspk;
                    stats.trialtype.(trialtypes{i})(j).events.reward.time = ts;
                    stats.trialtype.(trialtypes{i})(j).events.reward.peakresp = ...         % significance of peak response
                        EvaluatePeakresponse(trials_spks_temp2,prs.ts.reward,temporal_binwidth,peaktimewindow,minpeakprominence,nbootstraps,mintrialsforstats);
                end
                %% aligned to saccades
                if any(strcmp(gettuning,'saccade'))
                    trials_spks_temp2 = ShiftSpikes4saccades(trials_spks_temp,{events_temp.t_sac});
                    [nspk,ts] = Spiketimes2Rate(trials_spks_temp2,prs.ts.saccade,temporal_binwidth);
                    stats.trialtype.(trialtypes{i})(j).events.saccade.rate = nspk;
                    stats.trialtype.(trialtypes{i})(j).events.saccade.time = ts;
                    stats.trialtype.(trialtypes{i})(j).events.saccade.peakresp = ...         % significance of peak response
                        EvaluatePeakresponse(trials_spks_temp2,prs.ts.reward,temporal_binwidth,peaktimewindow,minpeakprominence,nbootstraps,mintrialsforstats);
                end
            end
        end
    end
end

%% cross-correlation and tuning to continuous variables (requires nonparametric-regression package: https://github.com/kaushik-l/nonparametric-regression.git)
if compute_tuning
    gettuning = [prs.tuning_continuous {'r_stop','eye_verhor','sacmag','sacdir','sacmagdir','targ_ver','targ_hor','targ_verhor','tte_ver','tte_hor','tte_mag','tte_verhor','r','theta'}];
    for i=1:length(trialtypes) % compute tuning curves using all trials, rather than separately for each condition
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        fprintf(['.........estimating tuning curves :: trialtype: ' (trialtypes{i}) '\n']);
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).continuous = stats.trialtype.all.continuous;
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                trials_spks_temp = trials_spks(trlindx);
                %% define time windows for computing tuning
                timewindow_move = [[events_temp.t_move]' [events_temp.t_stop]']; % when the subject is moving
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                timewindow_Toff = [[events_temp.t_targ]'+0.3 [events_temp.t_stop]']; % when the target has disappeared
                %% linear velocity, v
                if any(strcmp(gettuning,'v'))
                    stats.trialtype.(trialtypes{i})(j).continuous.v = ...
                        ComputeTuning({continuous_temp.v},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.v);
                end
                %% angular velocity, w
                if any(strcmp(gettuning,'w'))
                    stats.trialtype.(trialtypes{i})(j).continuous.w = ...
                        ComputeTuning({continuous_temp.w},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.w);
                end
                %% velocity, vw (two dimensional)
                if any(strcmp(gettuning,'vw'))
                    stats.trialtype.(trialtypes{i})(j).continuous.vw = ...
                        ComputeTuning2D({continuous_temp.v},{continuous_temp.w},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,prs.tuning,prs.tuning_method);
                end
                %% linear acceleration, a
                if any(strcmp(gettuning,'r_accel'))
                    a = cellfun(@(x) diff(x)/dt,{continuous_temp.v},'UniformOutput',false);
                    a_ts = cellfun(@(x) x(2:end),{continuous_temp.ts},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.a = ...
                        ComputeTuning(a,a_ts,{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,[]);
                end
                %% angular acceleration, alpha
                if any(strcmp(gettuning,'theta_accel'))
                    alpha = cellfun(@(x) diff(x)/dt,{continuous_temp.w},'UniformOutput',false);
                    alpha_ts = cellfun(@(x) x(2:end),{continuous_temp.ts},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.alpha = ...
                        ComputeTuning(alpha,alpha_ts,{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,[]);
                end
                %% acceleration, aalpha (two dimensional)
                if any(strcmp(gettuning,'rtheta_accel'))
                    stats.trialtype.(trialtypes{i})(j).continuous.aalpha = ...
                        ComputeTuning2D(a,alpha,a_ts,{trials_spks_temp.tspk},timewindow_move,prs.tuning,prs.tuning_method);
                end
                %% magnitude of linear velocity, |v|
                if any(strcmp(gettuning,'v_abs'))
                    v_abs = cellfun(@abs,{continuous_temp.v},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.v_abs = ...
                        ComputeTuning(v_abs,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.v);
                end
                %% magnitude of angular velocity, |w|
                if any(strcmp(gettuning,'w_abs'))
                    w_abs = cellfun(@abs,{continuous_temp.w},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.w_abs = ...
                        ComputeTuning(w_abs,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,[]);
                end
                %% magnitude of linear acceleration, |a|
                if any(strcmp(gettuning,'a_abs'))
                    a_abs = cellfun(@(x) abs(diff(x)/dt),{continuous_temp.v},'UniformOutput',false);
                    a_abs_ts = cellfun(@(x) x(2:end),{continuous_temp.ts},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.a_abs = ...
                        ComputeTuning(a_abs,a_abs_ts,{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,[]);
                end
                %% magnitude of angular acceleration, |alpha|
                if any(strcmp(gettuning,'alpha_abs'))
                    alpha_abs = cellfun(@(x) abs(diff(x)/dt),{continuous_temp.w},'UniformOutput',false);
                    alpha_abs_ts = cellfun(@(x) x(2:end),{continuous_temp.ts},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.alpha_abs = ...
                        ComputeTuning(alpha_abs,alpha_abs_ts,{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,[]);
                end
                %% vertical eye position, veye
                if any(strcmp(gettuning,'eye_ver'))
                    veye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false); % average both eyes (if available)
                    veye = cellfun(@(x) x(:), veye,'un',0); 
                    for m = 1:numel(veye); if sum(veye{m}==0)==numel(veye{m}); veye{m}(:) = nan; end; end
                    stats.trialtype.(trialtypes{i})(j).continuous.veye = ...
                        ComputeTuning(veye,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.eye_ver);
                end
                %% horizontal eye position, heye
                if any(strcmp(gettuning,'eye_hor'))
                    heye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false); % average both eyes (if available)
                    heye = cellfun(@(x) x(:), heye,'un',0);
                    for m = 1:numel(heye); if sum(heye{m}==0)==numel(heye{m}); heye{m}(:) = nan; end; end
                    stats.trialtype.(trialtypes{i})(j).continuous.heye = ...
                        ComputeTuning(heye,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.eye_hor);
                end
                %% eye position, eye_verhor (two dimensional)
                if any(strcmp(gettuning,'eye_verhor'))
                    stats.trialtype.(trialtypes{i})(j).continuous.vheye = ...
                        ComputeTuning2D(veye,heye,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_move,prs.tuning,prs.tuning_method);
                end
                %% vertical target position on screen, targ_ver
                if any(strcmp(gettuning,'targ_ver'))
                    vtarg = cellfun(@(x) x(:)',{continuous_temp.z_targ},'UniformOutput',false); % average both eyes (if available)
                    vtarg = cellfun(@(x) x(:), vtarg,'un',0);
                    stats.trialtype.(trialtypes{i})(j).continuous.vtarg = ...
                        ComputeTuning(vtarg,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.eye_ver);
                end
                %% horizontal target position on screen, targ_ver
                if any(strcmp(gettuning,'targ_hor'))
                    htarg = cellfun(@(x) x(:)',{continuous_temp.z_targ},'UniformOutput',false); % average both eyes (if available)
                    htarg = cellfun(@(x) x(:), htarg,'un',0);
                    stats.trialtype.(trialtypes{i})(j).continuous.htarg = ...
                        ComputeTuning(htarg,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.eye_ver);
                end
                %% target position, targ_verhor (two dimensional)
                if any(strcmp(gettuning,'targ_verhor'))
                    stats.trialtype.(trialtypes{i})(j).continuous.vhtarg = ...
                        ComputeTuning2D(vtarg,htarg,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,prs.tuning,prs.tuning_method);
                end
                %% vertical target-tracking error, tte_ver
                if any(strcmp(gettuning,'tte_ver'))
                    vtte = cellfun(@(x,y) x - y, vtarg,veye,'UniformOutput',false); 
                    stats.trialtype.(trialtypes{i})(j).continuous.vtte = ...
                        ComputeTuning(vtte,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.eye_hor);
                end
                %% horizontal target-tracking error, tte_hor
                if any(strcmp(gettuning,'tte_hor'))
                    htte = cellfun(@(x,y) x - y, htarg,heye,'UniformOutput',false); 
                    stats.trialtype.(trialtypes{i})(j).continuous.htte = ...
                        ComputeTuning(htte,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.eye_hor);
                end
                %% magnitude of target-tracking error, tte_mag
                if any(strcmp(gettuning,'tte_mag'))
                    tte_mag = cellfun(@(x,y) sqrt(x.^2 + y.^2),htte,vtte,'un',0);
                    stats.trialtype.(trialtypes{i})(j).continuous.ttemag = ...
                        ComputeTuning(tte_mag,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.sac_mag);
                end
                %% target-tracking error, tte_verhor (two dimensional)
                if any(strcmp(gettuning,'tte_verhor'))
                    stats.trialtype.(trialtypes{i})(j).continuous.vhtte = ...
                        ComputeTuning2D(vtte,htte,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_Toff,prs.tuning,prs.tuning_method);
                end
                %% displacement, r
                if any(strcmp(gettuning,'r'))
                    r = cellfun(@(x,y) sqrt((x(:)-x0).^2 + (y(:)-y0).^2),{continuous_temp.xmp},{continuous_temp.ymp},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.r = ...
                        ComputeTuning(r,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.r_targ);
                end
                %% bearing, theta
                if any(strcmp(gettuning,'theta'))
                    theta = cellfun(@(x,y) atan2d(x(:)-x0,y(:)-y0),{continuous_temp.xmp},{continuous_temp.ymp},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.theta = ...
                        ComputeTuning(theta,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.theta_targ);
                end
                %% position, rtheta (two dimensional)
                if any(strcmp(gettuning,'rtheta'))
                    stats.trialtype.(trialtypes{i})(j).continuous.rtheta = ...
                        ComputeTuning2D(r,theta,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,prs.tuning,prs.tuning_method);
                end
                %% distance, d (refine -- use t_targ instead of 0?)
                if any(strcmp(gettuning,'d'))
                    d = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.d = ...
                        ComputeTuning(d,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.d);
                end
                %% heading, phi
                if any(strcmp(gettuning,'phi'))
                    phi = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                    stats.trialtype.(trialtypes{i})(j).continuous.phi = ...
                        ComputeTuning(phi,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.phi);
                end
                %% path, dphi (two dimensional)
                if any(strcmp(gettuning,'dphi'))
                    stats.trialtype.(trialtypes{i})(j).continuous.dphi = ...
                        ComputeTuning2D(d,phi,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,prs.tuning,prs.tuning_method);
                end
                %% distance to target, r_targ
                if any(strcmp(gettuning,'r_targ'))
                    r_targ = behv_stats.pos_rel.r_targ(trlindx);
                    stats.trialtype.(trialtypes{i})(j).continuous.r_targ = ...
                        ComputeTuning(r_targ,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.r_targ);
                end
                %% angle to target, theta_targ
                if any(strcmp(gettuning,'theta_targ'))
                    theta_targ = behv_stats.pos_rel.theta_targ(trlindx);
                    stats.trialtype.(trialtypes{i})(j).continuous.theta_targ = ...
                        ComputeTuning(theta_targ,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.theta_targ);
                end
                %% distance to stop, r_stop
                if any(strcmp(gettuning,'r_stop'))
                    r_stop = behv_stats.pos_rel.r_stop(trlindx);
                    stats.trialtype.(trialtypes{i})(j).continuous.r_stop = ...
                        ComputeTuning(r_stop,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.r_targ);
                end
                %% angle to stop, theta_stop
                if any(strcmp(gettuning,'r_stop'))
                    theta_stop = behv_stats.pos_rel.theta_stop(trlindx);
                    stats.trialtype.(trialtypes{i})(j).continuous.theta_stop = ...
                        ComputeTuning(theta_stop,{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.theta_targ);
                end
                %% saccade direction
                if any(strcmp(gettuning,'sacdir'))
                    sac_dir = arrayfun(@(x) {x.sac_dir},events_temp);
                    sac_t = arrayfun(@(x) {x.t_sac},events_temp);
                    perievent_t = prs.ts.saccade(end)-prs.ts.saccade(1);
                    trials_spks_temp2 = ShiftSpikes4saccades(trials_spks_temp,{events_temp.t_sac});
                    t_spks = arrayfun(@(x) {x.tspk >= prs.ts.saccade(1) & x.tspk <= prs.ts.saccade(end)},trials_spks_temp2);
                                        
                    stats.trialtype.(trialtypes{i})(j).continuous.sac_dir = ...
                        ComputeTuning_Incontinuous(sac_dir,sac_t,t_spks,timewindow_path,perievent_t,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.sac_dir);
                end
                %% saccade magnitude
                if any(strcmp(gettuning,'sacmag'))
                    sac_mag = arrayfun(@(x) {x.sac_mag},events_temp);
                    stats.trialtype.(trialtypes{i})(j).continuous.sac_mag = ...
                        ComputeTuning_Incontinuous(sac_mag,sac_t,t_spks,timewindow_path,perievent_t,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.sac_mag);
                end
                %% saccade direction and magnitude (two dimensional)
                if any(strcmp(gettuning,'sacmagdir'))
                    stats.trialtype.(trialtypes{i})(j).continuous.sac_magdir = ...
                        ComputeTuning2D_Incontinuous(sac_mag,sac_dir,sac_t,t_spks,timewindow_move,prs.tuning,prs.tuning_method);
                end
            end
        end
    end
end

%% fit generalised additive model to determine tuning to task variables (requires neuroGAM package: https://github.com/kaushik-l/neuroGAM.git)
if fitGAM_tuning
    GAM_prs.varname = prs.GAM_varname; varname = GAM_prs.varname;
    GAM_prs.vartype = prs.GAM_vartype; vartype = GAM_prs.vartype;
    GAM_prs.basistype = prs.GAM_basistype; % enable only for branch smoothbasis of BuildGAM
    GAM_prs.nbins = prs.GAM_nbins;
    GAM_prs.binrange = [];
    GAM_prs.nfolds = prs.nfolds;
    GAM_prs.dt = prs.dt;
    GAM_prs.filtwidth = prs.neuralfiltwidth;
    GAM_prs.linkfunc = prs.GAM_linkfunc;
    GAM_prs.lambda = prs.GAM_lambda;
    GAM_prs.alpha = prs.GAM_alpha;
    GAM_prs.varchoose = prs.GAM_varchoose;
    GAM_prs.method = prs.GAM_method;
    for i=1% if i=1, fit model using data from all trials rather than separately to data from each condition
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        fprintf(['.........fitting GAM model :: trialtype: ' (trialtypes{i}) '\n']);
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).GAM.(GAM_prs.linkfunc) = stats.trialtype.all.GAM.(GAM_prs.linkfunc);
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                trials_spks_temp = trials_spks(trlindx);
                if ~isempty(lfps), trials_lfps_temp = lfps(prs.channel_id).trials(trlindx); end
                %% select variables of interest and load their details
                vars = cell(length(varname),1);
                GAM_prs.binrange = cell(1,length(varname));
                for k=1:length(varname)
                    if isfield(continuous_temp,varname(k)), vars{k} = {continuous_temp.(varname{k})};
                    elseif isfield(behv_stats.pos_rel,varname(k)), vars{k} = behv_stats.pos_rel.(varname{k})(trlindx);
                    elseif strcmp(varname(k),'d')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'phi')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'eye_ver')
                        isnan_le = all(isnan(cell2mat({continuous_temp.zle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.zre}')));
                        if isnan_le, vars{k} = {continuous_temp.zre};
                        elseif isnan_re, vars{k} = {continuous_temp.zle};
                        else vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'eye_hor')
                        isnan_le = all(isnan(cell2mat({continuous_temp.yle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.yre}')));
                        if isnan_le, vars{k} = {continuous_temp.yre};
                        elseif isnan_re, vars{k} = {continuous_temp.yle};
                        else vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'phase')
                        vars{k} = cellfun(@(x) angle(hilbert(x)), {trials_lfps_temp.lfp},'UniformOutput',false);
                    elseif strcmp(vartype(k),'event')
                        if ~strcmp(varname(k),'spikehist'), vars{k} = [events_temp.(prs.varlookup(varname{k}))]; else, vars{k} = []; end
                        if strcmp(varname(k),'target_OFF'), vars{k} = vars{k} + prs.fly_ONduration; end % target_OFF = t_targ + fly_ONduration
                    end
                    GAM_prs.binrange{k} = prs.binrange.(varname{k});
                end
                %% define time windows for computing tuning
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                timewindow_full = [min([events_temp.t_move],[events_temp.t_targ]) - prs.pretrial ;... % from "min(move,targ) - pretrial_buffer"
                    [events_temp.t_end] + prs.posttrial]'; % till "end + posttrial_buffer"
                %% concatenate data from all trials
                xt = []; yt = [];
                for k=1:length(vars)
                    if any(strcmp(varname(k),{'h1','h2'}))
                        [xt(:,k),~,yt] = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                        xt(:,k) = xt(:,k)/100; % 1 unit = 100pixels/s
                    elseif ~strcmp(vartype(k),'event')
                        [xt(:,k),~,yt] = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    elseif ~strcmp(varname(k),'spikehist')                        
                        [~,xt(:,k),yt] = ConcatenateTrials([],mat2cell(vars{k}',ones(length(events_temp),1)),{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    end
                end
                if any(strcmp(varname,'spikehist')), xt(:,strcmp(varname,'spikehist')) = yt; end % pass spike train back as an input to fit spike-history kernel
                %% model fitting and selection
                xt = mat2cell(xt,size(xt,1),ones(1,size(xt,2))); % convert to cell
                models = BuildGAM(xt,yt,GAM_prs);
                stats.trialtype.(trialtypes{i})(j).GAM.(GAM_prs.linkfunc) = models;
            end
        end
    end
    %% estimate variance explained by GAM
    if GAM_varexp
        if ~all(GAM_prs.varchoose == 0), warning('cannot compute variance explained if varchoose = 1');
        elseif ~strcmpi(GAM_prs.method,'fastbackward'), warning('set "fastbackward" as your model selection method for computing variance explained');
        else
            nvars = numel(GAM_prs.varname);
            modelclasses = cell2mat(stats.trialtype.all.GAM.log.class);            
            varexp_full = models.trainFit{all(modelclasses == ones(1,nvars),2)}(:,1);
%             if any(~strcmp(vartype,'event'))
%                 GAM_prs.method = 'FastForward'; models_fwd = BuildGAM(xt,yt,GAM_prs); 
%                 modelclasses_fwd = cell2mat(models_fwd.class); 
%             end
            for i=1:nvars
                % for var explained by discrete events, use the reduction in varexp when removing that variable
%                 if strcmpi(vartype{i},'event')
                    indx = true(1,nvars); indx(i) = false;
                    stats.trialtype.all.GAM.log.varexp{i} = mean(varexp_full - models.trainFit{all(modelclasses == indx,2)}(:,1));
                    if strcmpi(vartype{i},'event')
                        t_frac = (sum(~isnan(vars{i}))*diff(GAM_prs.binrange{i}))/(dt*numel(xt{i}));
                        stats.trialtype.all.GAM.log.varexp{i} = stats.trialtype.all.GAM.log.varexp{i}/t_frac;
                    end
                % for var explained by continuous vars, use the varexp with just that variable as predictor
%                 elseif any(strcmpi(vartype{i},{'1D','2D'}))
%                     indx = false(1,nvars); indx(i) = true;
%                     stats.trialtype.all.GAM.log.varexp{i} = mean(models_fwd.trainFit{all(modelclasses_fwd == indx,2)}(:,1));
%                 end
            end
        end
    end
end

%% Fit a 2-layer neural network model (requires spykesML Python package: https://github.com/KordingLab/spykesML.git)
if fitNNM
    NNM_prs.varname = prs.NNM_varname; varname = NNM_prs.varname;
    NNM_prs.vartype = prs.NNM_vartype; vartype = NNM_prs.vartype;
    NNM_prs.nbins = prs.NNM_nbins;
    NNM_prs.binrange = [];
    NNM_prs.filtwidth = prs.neuralfiltwidth;
    NNM_prs.method = prs.NNM_method;
    nvars = numel(NNM_prs.varname);
    for i=1% if i=1, fit model using data from all trials rather than separately to data from each condition
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        fprintf(['.........fitting Neural network model :: trialtype: ' (trialtypes{i}) '\n']);
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).NNM.(NNM_prs.method) = stats.trialtype.all.NNM.(NNM_prs.method);
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                trials_spks_temp = trials_spks(trlindx);
                if ~isempty(lfps), trials_lfps_temp = lfps(prs.channel_id).trials(trlindx); end
                %% select variables of interest and load their details
                vars = cell(length(varname),1);
                NNM_prs.binrange = cell(1,length(varname));
                for k=1:length(varname)
                    if isfield(continuous_temp,varname(k)), vars{k} = {continuous_temp.(varname{k})};
                    elseif isfield(behv_stats.pos_rel,varname(k)), vars{k} = behv_stats.pos_rel.(varname{k})(trlindx);
                    elseif strcmp(varname(k),'d')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'phi')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'eye_ver')
                        isnan_le = all(isnan(cell2mat({continuous_temp.zle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.zre}')));
                        if isnan_le, vars{k} = {continuous_temp.zre};
                        elseif isnan_re, vars{k} = {continuous_temp.zle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'eye_hor')
                        isnan_le = all(isnan(cell2mat({continuous_temp.yle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.yre}')));
                        if isnan_le, vars{k} = {continuous_temp.yre};
                        elseif isnan_re, vars{k} = {continuous_temp.yle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'phase')
                        vars{k} = cellfun(@(x) angle(hilbert(x)), {trials_lfps_temp.lfp},'UniformOutput',false);
                    elseif strcmp(vartype(k),'event')
                        if ~strcmp(varname(k),'spikehist'), vars{k} = [events_temp.(prs.varlookup(varname{k}))]; else, vars{k} = []; end
                        if strcmp(varname(k),'target_OFF'), vars{k} = vars{k} + prs.fly_ONduration; end % target_OFF = t_targ + fly_ONduration
                    end
                    NNM_prs.binrange{k} = prs.binrange.(varname{k});
                end
                %% define time windows for computing tuning
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                timewindow_full = [min([events_temp.t_move],[events_temp.t_targ]) - prs.pretrial ;... % from "min(move,targ) - pretrial_buffer"
                    [events_temp.t_end] + prs.posttrial]'; % till "end + posttrial_buffer"
                %% concatenate data from all trials
                xt = []; yt = [];
                for k=1:length(vars)
                    if ~strcmp(vartype(k),'event')
                        [xt(:,k),~,yt] = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    elseif ~strcmp(varname(k),'spikehist')
                        [~,xt(:,k),yt] = ConcatenateTrials([],mat2cell(vars{k}',ones(length(events_temp),1)),{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    end
                end
                if any(strcmp(varname,'spikehist')), xt(:,strcmp(varname,'spikehist')) = yt; end % pass spike train back as an input to fit spike-history kernel
                %% model fitting and selection
                xt = mat2cell(xt,size(xt,1),ones(1,size(xt,2))); % convert to cell
                indx = find(strcmp(NNM_prs.vartype,'event')); for k=indx, NNM_prs.binrange{k} = round(NNM_prs.binrange{k}/dt); end
                for k=1:nvars, x{k} = Encode1hot(xt{k}, NNM_prs.vartype{k}, NNM_prs.binrange{k}, NNM_prs.nbins{k}); end % use 1-hot encoding for neurons in the input layer
                X = cell2mat(x); y = yt'; %(conv(yt,h,'same')')/dt;
                save('tempdata_Xy.mat','X','y');  % store data in a .mat file for scipy to access
                [messenger,model] = system(['python ' 'C:\Users\jklakshm\Documents\GitHub\firefly-monkey\Analyse\FitNNmodel.py']); % switch to command line to run Python code
                delete('tempdata_Xy.mat'); % destroy the .mat file
                if ~messenger, stats.trialtype.(trialtypes{i})(j).NNM.(NNM_prs.method) = model;
                else, warning(['python script failed with the following error: ' messenger]); end
            end
        end
    end
end

%% Spike-LFP analysis
if ~isempty(lfps)
if analyse_spikeLFPrelation
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    trials_lfps = lfps(prs.channel_id).trials; % only use LFP from the same channel on which the unit was recorded
    % phase analysis
    for i=1:ntrls, trials_lfps(i).phase = angle(hilbert(trials_lfps(i).lfp)); end
    for i=1:length(trialtypes) % i=1:2  => 'all' and 'reward'
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).continuous.lfps = stats.trialtype.all.continuous.lfps;
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                trials_spks_temp = trials_spks(trlindx);
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                trials_lfps_temp = trials_lfps(trlindx);
                %% define time windows for computing tuning
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                %% compute tuning to phase
                stats.trialtype.(trialtypes{i})(j).continuous.lfps.phase = ...
                    ComputeTuning({trials_lfps_temp.phase},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,[],prs.tuning,prs.tuning_method);
                %% compute tuning to phase and v
%                 stats.trialtype.(trialtypes{i})(j).continuous.lfps.vphase = ...
%                     ComputeTuning2D({continuous_temp.v},{trials_lfps_temp.phase},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,prs.tuning,prs.tuning_method);
                %% compute tuning to phase and w
%                 stats.trialtype.(trialtypes{i})(j).continuous.lfps.wphase = ...
%                     ComputeTuning2D({continuous_temp.w},{trials_lfps_temp.phase},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,prs.tuning,prs.tuning_method);
                %% spike-triggered average of LFP
                stats.trialtype.(trialtypes{i})(j).continuous.lfps.sta = SpikeTriggeredLFP({trials_lfps_temp.lfp},{continuous_temp.ts},{trials_spks_temp.tspk},...
                    timewindow_path,sta_window,duration_nanpad,spectralparams);
            end
        end
    end
end

if analyse_temporalphase
    trials_lfps = lfps(prs.channel_id).trials; % only use LFP from the same channel on which the unit was recorded
    for i=1:ntrls, trials_lfps(i).phase = angle(hilbert(trials_lfps(i).lfp)); end
    trlindx = behv_stats.trialtype.all.trlindx;
    trials_spks_temp = trials_spks(trlindx);
    events_temp = events(trlindx);
    continuous_temp = continuous(trlindx);
    trials_lfps_temp = trials_lfps(trlindx);
    t_stop =  [events_temp.t_stop]'; % when the subject is moving
    prs_temp = prs; prs_temp.tuning.nbins1d_binning = prs.num_phasebins; prs_temp.tuning_method = 'local-linear'; % use different binning just for this analysis
    for l = 1:length(phase_slidingwindow)-1
        %% define time windows for computing tuning
        timewindow_use = [[events_temp.t_targ]'+phase_slidingwindow(l) [events_temp.t_targ]'+phase_slidingwindow(l+1)]; % window to use for computing phase
        trlindx2 = timewindow_use(:,2) < t_stop;
        %% compute tuning to phase        
        stats.trialtype.all.continuous.temporalphase.phi(l) = ...
            ComputeTuning({trials_lfps_temp(trlindx2).phase},{continuous_temp(trlindx2).ts},{trials_spks_temp(trlindx2).tspk},timewindow_use(trlindx2,:),duration_zeropad,corr_lag,[],prs_temp.tuning,prs_temp.tuning_method);
    end
    stats.trialtype.all.continuous.temporalphase.t = 0.5*(phase_slidingwindow(1:end-1) + phase_slidingwindow(2:end));
end

if analyse_spikeLFPrelation_allLFPs
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    nlfps = length(lfps);
    for k=1:nlfps
        fprintf(['                .....LFP ' num2str(k) '\n']);
        trials_lfps = lfps(k).trials;
        % phase analysis
        for i=1:ntrls
            trials_lfps(i).phase = angle(hilbert(trials_lfps(i).lfp));
        end
        for i=1:2%length(trialtypes) % i=1:2  => 'all' and 'reward'
            nconds = length(behv_stats.trialtype.(trialtypes{i}));
            if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
            for j=1:nconds
                if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                    stats.trialtype.(trialtypes{i})(j).continuous.lfps(k) = stats.trialtype.all.continuous.lfps(k);
                else
                    trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                    trials_spks_temp = trials_spks(trlindx);
                    events_temp = events(trlindx);
                    continuous_temp = continuous(trlindx);
                    trials_lfps_temp = trials_lfps(trlindx);
                    %% define time windows for computing tuning
                    timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                    %% compute tuning to phase
                    stats.trialtype.(trialtypes{i})(j).continuous.lfps(k).phase = ...
                        ComputeTuning({trials_lfps_temp.phase},{continuous_temp.ts},{trials_spks_temp.tspk},timewindow_path,duration_zeropad,corr_lag,[],prs.tuning,prs.tuning_method);
                    %% spike-triggered average of LFP
                    stats.trialtype.(trialtypes{i})(j).continuous.lfps(k).sta = SpikeTriggeredLFP({trials_lfps_temp.lfp},{continuous_temp.ts},{trials_spks_temp.tspk},...
                        timewindow_path,sta_window,duration_nanpad,spectralparams);
                end
            end
        end
    end
end

end
%% time-rescaling index

trlindx = behv_stats.trialtype.all.trlindx;
events_temp = events(trlindx);
trials_spks_temp = trials_spks(trlindx);
[stats.trialtype.all.intrinsic.scalingindex,stats.trialtype.all.intrinsic.lockingindex] = ...
    ComputeScalingindex(trials_spks_temp,events_temp,prs.ts_shortesttrialgroup,temporal_binwidth,prs.ntrialgroups);

