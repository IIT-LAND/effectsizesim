function heterogeneity_simulation(fname2save)
%
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

% plot population histograms
figure; set(gcf,'color','white');
subplot(3,1,1); hist(nonasd_data,100); grid on; xlim([-6 6]); title('Control population');
subplot(3,1,2); hist(pop_data_stacked,100); grid on; xlim([-6 6]); title('ASD population');
subplot(3,1,3); hist(pop_data_subgrp,500); grid on; xlim([-6 6]); title('ASD subtypes');
legend('ASD1','ASD2','ASD3','ASD4','ASD5');


%% Sample

% how many experiments to simulate?
n_exp = 10000;
effect_size20 = zeros(n_exp,1);
effect_size200 = zeros(n_exp,1);
effect_size2000 = zeros(n_exp,1);
subgrp_percentage20 = zeros(n_subgrp,n_exp);
subgrp_percentage200 = zeros(n_subgrp,n_exp);
subgrp_percentage2000 = zeros(n_subgrp,n_exp);

% run simulated experiments
parfor i = 1:n_exp
    disp(i);

    % what is the sample size?
    sample_size = 20;
    [asd_sample_data, idx] = datasample(pop_data_stacked, sample_size);
    nonasd_sample_data = datasample(nonasd_data, sample_size);
    effect_size20(i,1) = cohens_d(nonasd_sample_data, asd_sample_data,1,0);
    subgrp_count_table = zeros(n_subgrp,3);
    subgrp_count_table(:,1) = 1:n_subgrp;
    for isub = 1:n_subgrp
        tmp_lab = subgrp_labels_stacked(idx);
        subgrp_count_table(isub,2) = sum(tmp_lab==isub);
        subgrp_count_table(isub,3) = sum(tmp_lab==isub)./sample_size;
    end
    subgrp_percentage20(:,i) = subgrp_count_table(:,3);

    % what is the sample size?
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

    % what is the sample size?
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

end % for


figure; set(gcf,'color','white');
subplot(3,2,1);
hist(effect_size20,100); grid on; xlabel('Standardized Effect Size');
ylabel('Count');
xlim([-2 2]);
title('Case-Control Effect: n=20');

subplot(3,2,2);
hist(subgrp_percentage20'); grid on; xlabel('Percentage of Sample');
ylabel('Count');
legend('ASD1','ASD2','ASD3','ASD4','ASD5');
xlim([0 0.6]);
title('Subtype Sampled: n=20');

subplot(3,2,3); set(gcf,'color','white');
hist(effect_size200,100); grid on; xlabel('Standardized Effect Size');
ylabel('Count');
xlim([-2 2]);
title('Case-Control Effect: n=200');

subplot(3,2,4);
hist(subgrp_percentage200'); grid on; xlabel('Percentage of Sample');
ylabel('Count');
legend('ASD1','ASD2','ASD3','ASD4','ASD5');
xlim([0 0.6]);
title('Subtype Sampled: n=200');

subplot(3,2,5); set(gcf,'color','white');
hist(effect_size2000,100); grid on; xlabel('Standardized Effect Size');
ylabel('Count');
xlim([-2 2]);
title('Case-Control Effect: n=2000');

subplot(3,2,6);
hist(subgrp_percentage2000'); grid on; xlabel('Percentage of Sample');
ylabel('Count');
legend('ASD1','ASD2','ASD3','ASD4','ASD5');
xlim([0 0.6]);
title('Subtype Sampled: n=2000');




% plot population histograms
figure; set(gcf,'color','white');
subplot(3,2,1); hist(nonasd_data,100); grid on; xlim([-6 6]); title('Control population'); xlabel('DV'); ylabel('Count');
subplot(3,2,3); hist(pop_data_stacked,100); grid on; xlim([-6 6]); title('ASD population'); xlabel('DV'); ylabel('Count');
subplot(3,2,5); hist(pop_data_subgrp,200); grid on; xlim([-6 6]); title('ASD subtypes'); xlabel('DV'); ylabel('Count');
%legend('ASD1','ASD2','ASD3','ASD4','ASD5');

subplot(3,2,2);
hist(subgrp_percentage20'); grid on; xlabel('Sample Prevalence');
ylabel('Count');
% legend('ASD1','ASD2','ASD3','ASD4','ASD5');
xlim([0 0.6]);
title('n=20');

subplot(3,2,4);
hist(subgrp_percentage200'); grid on; xlabel('Sample Prevalence');
ylabel('Count');
% legend('ASD1','ASD2','ASD3','ASD4','ASD5');
xlim([0 0.6]);
title('n=200');

subplot(3,2,6);
hist(subgrp_percentage2000'); grid on; xlabel('Sample Prevalence');
ylabel('Count');
% legend('ASD1','ASD2','ASD3','ASD4','ASD5');
xlim([0 0.6]);
title('n=2000');


%% save plot
if ~isemtpy(fname2save)
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end % if ~isempty
