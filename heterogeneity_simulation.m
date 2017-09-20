function heterogeneity_simulation(fname2save)
%
%   Here we will simulate ASD population with 5 subtypes. Each subtype has a 20%
%   prevalence in the population. Both ASD and TD populations will be set to the
%   same mean and sd, so there will be no overall case-control difference. We
%   then simulate experiments with different sample sizes, and compute what is
%   the sample prevalence for the 5 ASD subtypes.
%
%   INPUT
%   fname2save = filename for pdf plot to save (leave empty if you don't want to save)
%


%% Population
% how many subgroups in the population
n_subgrp = 5;
% sample size of each subgroup
subgrp_pop_n = 200000;
% total population n
pop_n = n_subgrp*subgrp_pop_n;
% population mean for each subgroup
mu_subgrp = [-1, -0.5, 0, 0.5, 1];
% population standard deviation
pop_sd = 1;

% random normally distributed data for each subgroup with specific mean
pop_data_subgrp = [normrnd(mu_subgrp(1),pop_sd,subgrp_pop_n,1) ...
                   normrnd(mu_subgrp(2),pop_sd,subgrp_pop_n,1) ...
                   normrnd(mu_subgrp(3),pop_sd,subgrp_pop_n,1) ...
                   normrnd(mu_subgrp(4),pop_sd,subgrp_pop_n,1) ...
                   normrnd(mu_subgrp(5),pop_sd,subgrp_pop_n,1)];

% stack data into one vector
pop_data_stacked = reshape(pop_data_subgrp, pop_n, 1);
subgrp_labels_stacked = [ones(subgrp_pop_n,1); ones(subgrp_pop_n,1).*2; ...
                         ones(subgrp_pop_n,1).*3; ones(subgrp_pop_n,1).*4; ...
                         ones(subgrp_pop_n,1).*5];

% non-ASD comparison population with same mean
mean_asd = mean(pop_data_stacked);
% non-ASD SD is half that of ASD
nonasd_sd = std(pop_data_stacked)./2;
% generate nonasd_data
nonasd_data = normrnd(mean_asd, nonasd_sd, pop_n, 1);


%% Sample

% how many experiments to simulate?
n_exp = 10000;
% pre-allocate variables for the loop
% effect sizes
effect_size20 = zeros(n_exp,1);
effect_size200 = zeros(n_exp,1);
effect_size2000 = zeros(n_exp,1);
% sample prevalences
subgrp_percentage20 = zeros(n_subgrp,n_exp);
subgrp_percentage200 = zeros(n_subgrp,n_exp);
subgrp_percentage2000 = zeros(n_subgrp,n_exp);

% run simulated experiments
parfor i = 1:n_exp
    % disp(i);

    % n=20 ---------------------------------------------------------------------
    sample_size = 20;
    % randomly sample from both populations
    [asd_sample_data, idx] = datasample(pop_data_stacked, sample_size);
    nonasd_sample_data = datasample(nonasd_data, sample_size);
    % compute effect size
    effect_size20(i,1) = cohens_d(nonasd_sample_data, asd_sample_data,1,0);
    % count numbers of individuals in each subtype within the sample
    subgrp_count_table = zeros(n_subgrp,3);
    subgrp_count_table(:,1) = 1:n_subgrp;
    for isub = 1:n_subgrp
        tmp_lab = subgrp_labels_stacked(idx);
        subgrp_count_table(isub,2) = sum(tmp_lab==isub);
        subgrp_count_table(isub,3) = sum(tmp_lab==isub)./sample_size;
    end
    % calculate sample prevalence for each subtype
    subgrp_percentage20(:,i) = subgrp_count_table(:,3);

    % n=200 --------------------------------------------------------------------
    sample_size = 200;
    [asd_sample_data, idx] = datasample(pop_data_stacked, sample_size);
    nonasd_sample_data = datasample(nonasd_data, sample_size);
    effect_size200(i,1) = cohens_d(nonasd_sample_data, asd_sample_data,1,0);
    subgrp_count_table = zeros(n_subgrp,3);
    subgrp_count_table(:,1) = 1:n_subgrp;
    for isub = 1:n_subgrp
        tmp_lab = subgrp_labels_stacked(idx);
        subgrp_count_table(isub,2) = sum(tmp_lab==isub);
        subgrp_count_table(isub,3) = sum(tmp_lab==isub)./sample_size;
    end
    subgrp_percentage200(:,i) = subgrp_count_table(:,3);

    % n=2000 -------------------------------------------------------------------
    sample_size = 2000;
    [asd_sample_data, idx] = datasample(pop_data_stacked, sample_size);
    nonasd_sample_data = datasample(nonasd_data, sample_size);
    effect_size2000(i,1) = cohens_d(nonasd_sample_data, asd_sample_data,1,0);
    subgrp_count_table = zeros(n_subgrp,3);
    subgrp_count_table(:,1) = 1:n_subgrp;
    for isub = 1:n_subgrp
        tmp_lab = subgrp_labels_stacked(idx);
        subgrp_count_table(isub,2) = sum(tmp_lab==isub);
        subgrp_count_table(isub,3) = sum(tmp_lab==isub)./sample_size;
    end
    subgrp_percentage2000(:,i) = subgrp_count_table(:,3);

end % parfor



%% plot population histograms
XLIM1 = [-6, 6];
XLIM2 = [0, 0.6];
figure; set(gcf,'color','white');
% plot population distributions
subplot(3,2,1); hist(nonasd_data,100); grid on; xlim(XLIM1); title('Control population'); xlabel('DV'); ylabel('Count');
subplot(3,2,3); hist(pop_data_stacked,100); grid on; xlim(XLIM1); title('ASD population'); xlabel('DV'); ylabel('Count');
subplot(3,2,5); hist(pop_data_subgrp,200); grid on; xlim(XLIM1); title('ASD subtypes'); xlabel('DV'); ylabel('Count');

% plot sample prevalences
subplot(3,2,2);
hist(subgrp_percentage20'); grid on; xlabel('Sample Prevalence');
ylabel('Count');
xlim(XLIM2);
title('n=20');

subplot(3,2,4);
hist(subgrp_percentage200'); grid on; xlabel('Sample Prevalence');
ylabel('Count');
xlim(XLIM2);
title('n=200');

subplot(3,2,6);
hist(subgrp_percentage2000'); grid on; xlabel('Sample Prevalence');
ylabel('Count');
xlim(XLIM2);
title('n=2000');


%% save plot
if ~isemtpy(fname2save)
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end % if ~isempty
