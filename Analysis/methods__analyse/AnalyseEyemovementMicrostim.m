function [ts_eye_movement,eye_movement] = AnalyseEyemovementMicrostim(x_fly,y_fly,zle,yle,zre,yre,t_sac,t_stop,t_stim,ts,trlerrors,prs)

%% prs
delta = prs.interoculardist/2;
zt = -prs.height;
saccade_duration = prs.saccade_duration;
fly_ONduration = prs.fly_ONduration;
Nboots = prs.bootstrap_trl;
factor_downsample = 1; %prs.factor_downsample; % downsampling factor for storage
ntrls = length(x_fly);
pretrial = prs.pretrial;
posttrial = prs.posttrial;

%% sort trials by error
[~,errorindx] = sort(trlerrors);

%% eye position immediately after the first saccade following target onset
for i=1:ntrls
    % identify time of target fixation
    sacstart = []; sacend = []; sacampli = [];
    t_sac2 = t_sac{i};
    sac_indx = t_sac{i}>0 & t_sac{i}<2*fly_ONduration;
    if any(sac_indx)
        t_sacs = t_sac{i}(sac_indx);
        for j=1:length(t_sacs)
            sacstart(j) = find(ts{i}>(t_sacs(j)), 1);
            sacend(j) = find(ts{i}>(t_sacs(j) + saccade_duration), 1);
            sacampli(j) = nanmean([sum(abs(zle{i}(sacstart(j)) - zle{i}(sacend(j)))^2 + abs(yle{i}(sacstart(j)) - yle{i}(sacend(j)))^2) ...
                sum(abs(zre{i}(sacstart(j)) - zre{i}(sacend(j)))^2 + abs(yre{i}(sacstart(j)) - yre{i}(sacend(j)))^2)]);
        end
        t_fix(i) = t_sacs(sacampli == max(sacampli)) + saccade_duration/2;
    else, t_fix(i) = 0 + saccade_duration/2; 
    end % if no saccade detected, assume monkey was already fixating on target
    % remove saccade periods from eye position data
    sacstart = []; sacend = []; 
    for j=1:length(t_sac2)
        sacstart(j) = find(ts{i}>(t_sac2(j) - saccade_duration/2), 1);
        sacend(j) = find(ts{i}>(t_sac2(j) + saccade_duration/2), 1);
        xt{i}(sacstart(j):sacend(j)) = nan;  % fly x - position
        yle{i}(sacstart(j):sacend(j)) = nan; % left eye horizontal position
        yre{i}(sacstart(j):sacend(j)) = nan; % right eye horizontal position
        yt{i}(sacstart(j):sacend(j)) = nan;  % fly y - position
        zle{i}(sacstart(j):sacend(j)) = nan; % left eye vertical position
        zre{i}(sacstart(j):sacend(j)) = nan; % right eye vertical position
    end
    prestim = 0.5; poststim = 1;
    % select data around microstimulation
    xt{i} = x_fly{i}(ts{i}>(t_stim(i)-prestim) & ts{i}<(t_stim(i)+poststim)); yt{i} = y_fly{i}(ts{i}>(t_stim(i)-prestim) & ts{i}<(t_stim(i)+poststim));
    xt{i}(isnan(xt{i})) = xt{i}(find(~isnan(xt{i}),1)); yt{i}(isnan(yt{i})) = yt{i}(find(~isnan(yt{i}),1));
    yle{i} = yle{i}(ts{i}>(t_stim(i)-prestim) & ts{i}<(t_stim(i)+poststim)); yre{i} = yre{i}(ts{i}>(t_stim(i)-prestim) & ts{i}<(t_stim(i)+poststim));
    zle{i} = zle{i}(ts{i}>(t_stim(i)-prestim) & ts{i}<(t_stim(i)+poststim)); zre{i} = zre{i}(ts{i}>(t_stim(i)-prestim) & ts{i}<(t_stim(i)+poststim));
    if t_stim(i)<0.75 % disregard trials when stimulation is too early
        xt{i}(:) = nan; yt{i}(:) = nan;
        yle{i}(:) = nan; yle{i}(:) = nan; yre{i}(:) = nan; yre{i}(:) = nan;
        zle{i}(:) = nan; zle{i}(:) = nan; zre{i}(:) = nan; zre{i}(:) = nan;
    end
    
    % ground truth prediction for eye position (if the monkey really followed the target)
    yle_pred{i} = atan2d(xt{i} + delta, sqrt(yt{i}.^2 + zt^2));
    yre_pred{i} = atan2d(xt{i} - delta, sqrt(yt{i}.^2 + zt^2));
    zle_pred{i} = atan2d(zt , sqrt(yt{i}.^2 + (xt{i} + delta).^2));
    zre_pred{i} = atan2d(zt , sqrt(yt{i}.^2 + (xt{i} - delta).^2));
    ver_mean_pred{i} = nanmean([zle_pred{i} , zre_pred{i}],2); % mean vertical eye position (of the two eyes)
    hor_mean_pred{i} = nanmean([yle_pred{i} , yre_pred{i}],2); % mean horizontal eye position
    ver_diff_pred{i} = 0.5*(zle_pred{i} - zre_pred{i}); % 0.5*difference between vertical eye positions (of the two eyes)
    hor_diff_pred{i} = 0.5*(yle_pred{i} - yre_pred{i}); % 0.5*difference between horizontal eye positions
    % actual eye position
    ver_mean{i} = nanmean([zle{i} , zre{i}],2); % mean vertical eye position (of the two eyes)
    hor_mean{i} = nanmean([yle{i} , yre{i}],2); % mean horizontal eye position
    ver_diff{i} = 0.5*(zle{i} - zre{i}); % 0.5*difference between vertical eye positions (of the two eyes)
    hor_diff{i} = 0.5*(yle{i} - yre{i}); % 0.5*difference between horizontal eye positions
    % fly position
    rt{i} = sqrt(xt{i}.^2 + yt{i}.^2);
    thetat{i} = atan2d(xt{i},yt{i});
end
ts_eye_movement = linspace(-prestim,poststim,length(xt{1}));

%% regression
for i=1:ntrls
    eye_movement.flypos.r{i} = decimate(rt{i},factor_downsample);
    eye_movement.flypos.theta{i} = decimate(thetat{i},factor_downsample);
    eye_movement.eyepos.pred.ver_mean.val{i} = decimate(ver_mean_pred{i},factor_downsample);
    eye_movement.eyepos.pred.hor_mean.val{i} = decimate(hor_mean_pred{i},factor_downsample);
    eye_movement.eyepos.pred.ver_diff.val{i} = decimate(ver_diff_pred{i},factor_downsample);
    eye_movement.eyepos.pred.hor_diff.val{i} = decimate(hor_diff_pred{i},factor_downsample);
    eye_movement.eyepos.true.ver_mean.val{i} = decimate(ver_mean{i},factor_downsample);
    eye_movement.eyepos.true.hor_mean.val{i} = decimate(hor_mean{i},factor_downsample);
    eye_movement.eyepos.true.ver_diff.val{i} = decimate(ver_diff{i},factor_downsample);
    eye_movement.eyepos.true.hor_diff.val{i} = decimate(hor_diff{i},factor_downsample);
end

Nt = max(cellfun(@(x) length(x),rt)); % max number of timepoints
% temporal correlation between fly position & predicted eye position -
% aligned to stimulation
rt1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],rt,'UniformOutput',false));
thetat1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],thetat,'UniformOutput',false));

ver_mean_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean_pred,'UniformOutput',false));
hor_mean_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean_pred,'UniformOutput',false));
ver_diff_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_diff_pred,'UniformOutput',false));
hor_diff_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_diff_pred,'UniformOutput',false));
[eye_movement.eyepos.pred.ver_mean.rho.r.stimaligned,eye_movement.eyepos.pred.ver_mean.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',ver_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.hor_mean.rho.r.stimaligned,eye_movement.eyepos.pred.hor_mean.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',hor_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.ver_diff.rho.r.stimaligned,eye_movement.eyepos.pred.ver_diff.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',ver_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.hor_diff.rho.r.stimaligned,eye_movement.eyepos.pred.hor_diff.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',hor_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.ver_mean.rho.theta.stimaligned,eye_movement.eyepos.pred.ver_mean.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',ver_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.hor_mean.rho.theta.stimaligned,eye_movement.eyepos.pred.hor_mean.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',hor_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.ver_diff.rho.theta.stimaligned,eye_movement.eyepos.pred.ver_diff.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',ver_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred.hor_diff.rho.theta.stimaligned,eye_movement.eyepos.pred.hor_diff.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',hor_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);

% temporal correlation between fly position & true eye position -
% aligned to target fixation
ver_mean1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean,'UniformOutput',false));
hor_mean1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean,'UniformOutput',false));
ver_diff1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_diff,'UniformOutput',false));
hor_diff1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_diff,'UniformOutput',false));
[eye_movement.eyepos.true.ver_mean.rho.r.stimaligned,eye_movement.eyepos.true.ver_mean.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',ver_mean1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.hor_mean.rho.r.stimaligned,eye_movement.eyepos.true.hor_mean.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',hor_mean1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.ver_diff.rho.r.stimaligned,eye_movement.eyepos.true.ver_diff.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',ver_diff1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.hor_diff.rho.r.stimaligned,eye_movement.eyepos.true.hor_diff.pval.r.stimaligned] = arrayfun(@(i) corr(rt1(i,:)',hor_diff1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.ver_mean.rho.theta.stimaligned,eye_movement.eyepos.true.ver_mean.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',ver_mean1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.hor_mean.rho.theta.stimaligned,eye_movement.eyepos.true.hor_mean.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',hor_mean1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.ver_diff.rho.theta.stimaligned,eye_movement.eyepos.true.ver_diff.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',ver_diff1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.true.hor_diff.rho.theta.stimaligned,eye_movement.eyepos.true.hor_diff.pval.theta.stimaligned] = arrayfun(@(i) corr(thetat1(i,:)',hor_diff1(i,:)','Type','Spearman','rows','complete'), 1:Nt);

% timecourse of cosine similarity between predicted & true eye position -
% aligned to target fixation
[eye_movement.eyepos.pred_vs_true.ver_mean.rho.stimaligned,eye_movement.eyepos.pred_vs_true.ver_mean.pval.stimaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred_vs_true.hor_mean.rho.stimaligned,eye_movement.eyepos.pred_vs_true.hor_mean.pval.stimaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred_vs_true.ver_diff.rho.stimaligned,eye_movement.eyepos.pred_vs_true.ver_diff.pval.stimaligned] = arrayfun(@(i) corr(ver_diff1(i,:)',ver_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred_vs_true.hor_diff.rho.stimaligned,eye_movement.eyepos.pred_vs_true.hor_diff.pval.stimaligned] = arrayfun(@(i) corr(hor_diff1(i,:)',hor_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
pred = permute(cat(3,ver_mean_pred1 , hor_mean_pred1 , hor_diff_pred1),[3 1 2]);
true = permute(cat(3,ver_mean1 , hor_mean1 , hor_diff1),[3 1 2]);
cos_similarity = nan(Nboots,Nt); cos_similarity_shuffled = nan(Nboots,Nt);
for i=1:Nboots
    randtrls = randsample(ntrls,ntrls,1); randtrls2 = randsample(ntrls,ntrls,1);
    cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls)), 1:Nt);
    cos_similarity_shuffled(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls2)), 1:Nt);
end
eye_movement.eyepos.pred_vs_true.cos_similarity.mu.stimaligned = mean(cos_similarity);
eye_movement.eyepos.pred_vs_true.cos_similarity.sem.stimaligned = std(cos_similarity);
eye_movement.eyepos.pred_vs_true.cos_similarity_shuffled.mu.stimaligned = mean(cos_similarity_shuffled);
eye_movement.eyepos.pred_vs_true.cos_similarity_shuffled.sem.stimaligned = std(cos_similarity_shuffled);

% timecourse of cosine similarity between predicted & true eye position -
% trials grouped by error magnitude
ngroups = 5;
ntrls_per_group = (ntrls - mod(ntrls,ngroups))/ngroups;
errorindx = errorindx(1:ntrls_per_group*ngroups);
cos_similarity = nan(ngroups,Nt);
for i=1:ngroups
    trlgroup = errorindx(ntrls_per_group*(i-1) + 1:ntrls_per_group*i);
    cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,trlgroup),true(:,i,trlgroup)), 1:Nt);
end
eye_movement.eyepos.pred_vs_true.cos_similarity_vs_error.mu.stimaligned = cos_similarity;

% timecourse of centered cosine similarity between predicted & true eye position -
% aligned to target fixation
[eye_movement.eyepos.pred_vs_true.ver_mean.rho.stimaligned,eye_movement.eyepos.pred_vs_true.ver_mean.pval.stimaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred_vs_true.hor_mean.rho.stimaligned,eye_movement.eyepos.pred_vs_true.hor_mean.pval.stimaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred_vs_true.ver_diff.rho.stimaligned,eye_movement.eyepos.pred_vs_true.ver_diff.pval.stimaligned] = arrayfun(@(i) corr(ver_diff1(i,:)',ver_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[eye_movement.eyepos.pred_vs_true.hor_diff.rho.stimaligned,eye_movement.eyepos.pred_vs_true.hor_diff.pval.stimaligned] = arrayfun(@(i) corr(hor_diff1(i,:)',hor_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% regression - aligned to target onset
beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_mean.beta10.stimaligned,eye_movement.eyepos.pred_vs_true.ver_mean.beta00.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_mean.beta10.stimaligned,eye_movement.eyepos.pred_vs_true.hor_mean.beta00.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_diff1(i,:)',[ver_diff_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_diff.beta10.stimaligned,eye_movement.eyepos.pred_vs_true.ver_diff.beta00.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_diff1(i,:)',[hor_diff_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_diff.beta10.stimaligned,eye_movement.eyepos.pred_vs_true.hor_diff.beta00.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_mean.beta1.stimaligned,eye_movement.eyepos.pred_vs_true.ver_mean.beta0.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_mean.beta1.stimaligned,eye_movement.eyepos.pred_vs_true.hor_mean.beta0.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_diff1(i,:)',[ver_diff_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_diff.beta1.stimaligned,eye_movement.eyepos.pred_vs_true.ver_diff.beta0.stimaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_diff1(i,:)',[hor_diff_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_diff.beta1.stimaligned,eye_movement.eyepos.pred_vs_true.hor_diff.beta0.stimaligned] = deal(beta(1,:),beta(2,:));
pred = permute(cat(3,ver_mean_pred1 , hor_mean_pred1 , hor_diff_pred1),[3 1 2]); pred = (pred - repmat(nanmean(pred,3),[1 1 ntrls]))./repmat(nanstd(pred,[],3),[1 1 ntrls]);
true = permute(cat(3,ver_mean1 , hor_mean1 , hor_diff1),[3 1 2]); true = (true - repmat(nanmean(true,3),[1 1 ntrls]))./repmat(nanstd(true,[],3),[1 1 ntrls]);
cntr_cos_similarity = nan(Nboots,Nt);
for i=1:Nboots
    randtrls = randsample(ntrls,ntrls,1);
    cntr_cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls)), 1:Nt);
end
eye_movement.eyepos.pred_vs_true.cntr_cos_similarity.mu.stimaligned = mean(cntr_cos_similarity);
eye_movement.eyepos.pred_vs_true.cntr_cos_similarity.sem.stimaligned = std(cntr_cos_similarity);