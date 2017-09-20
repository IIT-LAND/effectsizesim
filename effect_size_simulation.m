function [Results] = effect_size_simulation(popES, fname2save)
%
%   Here we will simulate many experiments were there is a case-control
%   difference of some size (popES). We will simulate those experiments at many
%   different sample sizes. We will also calculate effect size inflation for
%   experiments that pass a nominal p<0.05 threshold.
%
%   INPUT
%   popES = set to the desired population effect size (e.g., 0.5)
%   fname2save = filename for pdf plot to save (leave empty if you don't want to save)
%
%   OUTPUT
%   Results = structure filled with average effect size inflation estimates
%

% set seed for reproducibility
rng(1);

%% Population

% population parameters
pop_mean1 = popES;
pop_sd1 = 1;
pop_mean2 = 0;
pop_sd2 = 1;

% population size
pop_size = 10000000;

% generate population data
pop1_data = normrnd(pop_mean1, pop_sd1, pop_size,1);
pop2_data = normrnd(pop_mean2, pop_sd2, pop_size,1);

D = cohens_d(pop1_data, pop2_data, 1, 0);

%% Simulate experiments with varying sample sizes and estimate effect size
% number of experiments to simulate
n_exp = 10000;

% sample sizes
sample_sizes = [20, 50, 100, 200, 1000, 2000];


% pre-allocate sample effect size estimate variable
d = zeros(n_exp,length(sample_sizes)); % matrix of effect sizes
t = zeros(size(d)); % vector of t-states
p = zeros(size(d)); % vector of p-values
reject_h0 = zeros(size(d)); % logical vector with 1 for rejected H0, otherwise 0

% Simulate n experiments
for i = 1:n_exp

    % n=20 ---------------------------------------------------------------------
    samp1_data = datasample(pop1_data,sample_sizes(1)); % randomly sample group1
    samp2_data = datasample(pop2_data,sample_sizes(1)); % randomly sample group2
    d(i,1) = cohens_d(samp1_data, samp2_data, 1, 0); % compute effect size
    [H,P,CI,STATS] = ttest2(samp1_data,samp2_data); % run hypothesis test
    t(i,1) = STATS.tstat; % grab t-stat
    p(i,1) = P; % grab p-value
    reject_h0(i,1) = H; % grab indicator of whether H0 is rejected

    % n=50 ---------------------------------------------------------------------
    samp1_data = datasample(pop1_data,sample_sizes(2));
    samp2_data = datasample(pop2_data,sample_sizes(2));
    d(i,2) = cohens_d(samp1_data, samp2_data, 1, 0);
    [H,P,CI,STATS] = ttest2(samp1_data,samp2_data);
    t(i,2) = STATS.tstat;
    p(i,2) = P;
    reject_h0(i,2) = H;

    % n=100 --------------------------------------------------------------------
    samp1_data = datasample(pop1_data,sample_sizes(3));
    samp2_data = datasample(pop2_data,sample_sizes(3));
    d(i,3) = cohens_d(samp1_data, samp2_data, 1, 0);
    [H,P,CI,STATS] = ttest2(samp1_data,samp2_data);
    t(i,3) = STATS.tstat;
    p(i,3) = P;
    reject_h0(i,3) = H;

    % n=200 --------------------------------------------------------------------
    samp1_data = datasample(pop1_data,sample_sizes(4));
    samp2_data = datasample(pop2_data,sample_sizes(4));
    d(i,4) = cohens_d(samp1_data, samp2_data, 1, 0);
    [H,P,CI,STATS] = ttest2(samp1_data,samp2_data);
    t(i,4) = STATS.tstat;
    p(i,4) = P;
    reject_h0(i,4) = H;

    % n=1000 -------------------------------------------------------------------
    samp1_data = datasample(pop1_data,sample_sizes(5));
    samp2_data = datasample(pop2_data,sample_sizes(5));
    d(i,5) = cohens_d(samp1_data, samp2_data, 1, 0);
    [H,P,CI,STATS] = ttest2(samp1_data,samp2_data);
    t(i,5) = STATS.tstat;
    p(i,5) = P;
    reject_h0(i,5) = H;

    % n=2000 -------------------------------------------------------------------
    samp1_data = datasample(pop1_data,sample_sizes(6));
    samp2_data = datasample(pop2_data,sample_sizes(6));
    d(i,6) = cohens_d(samp1_data, samp2_data, 1, 0);
    [H,P,CI,STATS] = ttest2(samp1_data,samp2_data);
    t(i,6) = STATS.tstat;
    p(i,6) = P;
    reject_h0(i,6) = H;

end % for i
reject_h0 = logical(reject_h0); % make reject_h0 a logical variable

%% figure out how inflated effect sizes are for rejected H0
for i = 1:size(d,2)
    % grab effect sizes from tests that had rejected H0
    data2use = d(reject_h0(:,i),i);

    % average effect size for rejected H0
    avg_inflate(i,1) = mean(data2use);

    % median effect size for rejected H0
    median_inflate(i,1) = median(data2use);

    % min effect size for rejected H0
    min_inflate(i,1) = min(data2use);

    % max effect size for rejected H0
    max_inflate(i,1) = max(data2use);

    % effect sizes for 99% of rejected H0
    ci_inflate(i,:) = prctile(data2use,[0.5 99.5]);

end % for i

%% calculate effect size inflation stats for rejected H0

% es_inflate_stats = [bottom 1%, mean, top 99%]
es_inflation_stats = [ci_inflate(:,1) avg_inflate ci_inflate(:,2)];
es_inflation_stats

% average effect size fold increase
avg_es_inflate_foldincrease = ((avg_inflate-D)./D);
% average effect size percentage increase
avg_es_inflate_percentincrease = avg_es_inflate_foldincrease.*100;
avg_es_inflate_percentincrease
% save avg_es_inflate_percentincrease to Results structure as output
Results.avg_es_inflate_percentincrease = avg_es_inflate_percentincrease;

% median effect size fold increase
median_es_inflate_foldincrease = ((median_inflate-D)./D);
% median effect size percentage increase
median_es_inflate_percentincrease = median_es_inflate_foldincrease.*100;
median_es_inflate_percentincrease

% max effect size fold increase
max_es_inflate_foldincrease = ((max_inflate-D)./D);
% average effect size percentage increase
max_es_inflate_percentincrease = max_es_inflate_foldincrease.*100;
max_es_inflate_percentincrease


%% Plot results
% number of bins for histogram
nbins = 100;
% xlimits to plot
XLIM = [-1 2];
lineWidth = 1; % how thick lines should be
pop_es_lineColor = [0 1 0]; % what color should population effect size line be
rejected_h0_color = [1 0 0]; % color for rejected_h0 distribution
faceAlpha = 1; % can modify transparency of rejected_H0 distribution here
% CI limits
ci_lim = [2.5 97.5];

% make figure, and set background to white
figure; set(gcf,'color','white');

% n=20 -------------------------------------------------------------------------
idx2use = 1;
subplot(3,2,idx2use);
% plot histogram of sample effect sizes
h = histogram(d(:,idx2use),nbins); h.FaceColor = ones(1,3)./2; h.EdgeColor = h.FaceColor; %axis square;
% plot histogram of effect sizes for rejected H0
hold on;
h = histogram(d(reject_h0(:,idx2use),idx2use),nbins); h.FaceColor = rejected_h0_color; h.EdgeColor = h.FaceColor; h.FaceAlpha = faceAlpha; %axis square;
% plot line for population effect
hold on; l = line([D D],ylim); l.Color = pop_es_lineColor; l.LineWidth = lineWidth;
% plot 95% confidence intervals for sample effect sizes
ci(idx2use,:) = prctile(d(:,idx2use),ci_lim);
l = line([ci(idx2use,1),ci(idx2use,1)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
l = line([ci(idx2use,2),ci(idx2use,2)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
grid on;
xlabel('Effect Size');
ylabel('Count');
xlim(XLIM);
title('n=20');

% n=50 -------------------------------------------------------------------------
idx2use = 2;
subplot(3,2,idx2use);
% plot histogram of sample effect sizes
h = histogram(d(:,idx2use),nbins); h.FaceColor = ones(1,3)./2; h.EdgeColor = h.FaceColor; %axis square;
% plot histogram of effect sizes for rejected H0
hold on;
h = histogram(d(reject_h0(:,idx2use),idx2use),nbins); h.FaceColor = rejected_h0_color; h.EdgeColor = h.FaceColor; h.FaceAlpha = faceAlpha; %axis square;
% plot line for population effect
hold on; l = line([D D],ylim); l.Color = pop_es_lineColor; l.LineWidth = lineWidth;
% plot 95% confidence intervals for sample effect sizes
ci(idx2use,:) = prctile(d(:,idx2use),ci_lim);
l = line([ci(idx2use,1),ci(idx2use,1)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
l = line([ci(idx2use,2),ci(idx2use,2)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
grid on;
xlabel('Effect Size');
ylabel('Count');
xlim(XLIM);
title('n=50');

% n=100 ------------------------------------------------------------------------
idx2use = 3;
subplot(3,2,idx2use);
% plot histogram of sample effect sizes
h = histogram(d(:,idx2use),nbins); h.FaceColor = ones(1,3)./2; h.EdgeColor = h.FaceColor; %axis square;
% plot histogram of effect sizes for rejected H0
hold on;
h = histogram(d(reject_h0(:,idx2use),idx2use),nbins); h.FaceColor = rejected_h0_color; h.EdgeColor = h.FaceColor; h.FaceAlpha = faceAlpha; %axis square;
% plot line for population effect
hold on; l = line([D D],ylim); l.Color = pop_es_lineColor; l.LineWidth = lineWidth;
% plot 95% confidence intervals for sample effect sizes
ci(idx2use,:) = prctile(d(:,idx2use),ci_lim);
l = line([ci(idx2use,1),ci(idx2use,1)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
l = line([ci(idx2use,2),ci(idx2use,2)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
grid on;
xlabel('Effect Size');
ylabel('Count');
xlim(XLIM);
title('n=100');

% n=200 ------------------------------------------------------------------------
idx2use = 4;
subplot(3,2,idx2use);
% plot histogram of sample effect sizes
h = histogram(d(:,idx2use),nbins); h.FaceColor = ones(1,3)./2; h.EdgeColor = h.FaceColor; %axis square;
% plot histogram of effect sizes for rejected H0
hold on;
h = histogram(d(reject_h0(:,idx2use),idx2use),nbins); h.FaceColor = rejected_h0_color; h.EdgeColor = h.FaceColor; h.FaceAlpha = faceAlpha; %axis square;
% plot line for population effect
hold on; l = line([D D],ylim); l.Color = pop_es_lineColor; l.LineWidth = lineWidth;
% plot 95% confidence intervals for sample effect sizes
ci(idx2use,:) = prctile(d(:,idx2use),ci_lim);
l = line([ci(idx2use,1),ci(idx2use,1)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
l = line([ci(idx2use,2),ci(idx2use,2)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
grid on;
xlabel('Effect Size');
ylabel('Count');
xlim(XLIM);
title('n=200');

% n=1000 -----------------------------------------------------------------------
idx2use = 5;
subplot(3,2,idx2use);
% plot histogram of sample effect sizes
h = histogram(d(:,idx2use),nbins); h.FaceColor = ones(1,3)./2; h.EdgeColor = h.FaceColor; %axis square;
% plot histogram of effect sizes for rejected H0
hold on;
h = histogram(d(reject_h0(:,idx2use),idx2use),nbins); h.FaceColor = rejected_h0_color; h.EdgeColor = h.FaceColor; h.FaceAlpha = faceAlpha; %axis square;
% plot line for population effect
hold on; l = line([D D],ylim); l.Color = pop_es_lineColor; l.LineWidth = lineWidth;
% plot 95% confidence intervals for sample effect sizes
ci(idx2use,:) = prctile(d(:,idx2use),ci_lim);
l = line([ci(idx2use,1),ci(idx2use,1)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
l = line([ci(idx2use,2),ci(idx2use,2)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
grid on;
xlabel('Effect Size');
ylabel('Count');
xlim(XLIM);
title('n=1000');

% n=2000 -----------------------------------------------------------------------
idx2use = 6;
subplot(3,2,idx2use);
% plot histogram of sample effect sizes
h = histogram(d(:,idx2use),nbins); h.FaceColor = ones(1,3)./2; h.EdgeColor = h.FaceColor; %axis square;
% plot histogram of effect sizes for rejected H0
hold on;
h = histogram(d(reject_h0(:,idx2use),idx2use),nbins); h.FaceColor = rejected_h0_color; h.EdgeColor = h.FaceColor; h.FaceAlpha = faceAlpha; %axis square;
% plot line for population effect
hold on; l = line([D D],ylim); l.Color = pop_es_lineColor; l.LineWidth = lineWidth;
% plot 95% confidence intervals for sample effect sizes
ci(idx2use,:) = prctile(d(:,idx2use),ci_lim);
l = line([ci(idx2use,1),ci(idx2use,1)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
l = line([ci(idx2use,2),ci(idx2use,2)],ylim); l.Color = zeros(1,3); l.LineWidth = lineWidth;
grid on;
xlabel('Effect Size');
ylabel('Count');
xlim(XLIM);
title('n=2000');

%% save plot
if ~isemtpy(fname2save)
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end % if ~isempty
