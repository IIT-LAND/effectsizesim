function effect_size_inflation_sim(fname2save)
%
%   Here we will calculate effect size inflation estimates over a range of
%   population effect sizes and sample sizes. We will then plot how effect size
%   inflation changes with different sample sizes and different population
%   effect sizes.
%
%   INPUT
%   fname2save = filename for pdf plot to save (leave empty if you don't want to save)
%

%% Simulate many experiments at different population effect sizes
% population effect sizes to loop over
d = 0:0.1:2;
% run a simulation at each effect size
for i = 1:length(d)
    Results = effect_size_simulation(d(i),'');
    % grab average effect size inflation estimates
    avg_es_inf(:,i) = Results.avg_es_inflate_percentincrease;
    close all; % close figure created by effect_size_simulation
end

%% make plot of effect size inflation at each sample size and population effect size
lineWidth = 2;

figure; set(gcf,'color','white');
p = plot(avg_es_inf(:,2:end)');
for i = 1:length(p); p(i).LineWidth = lineWidth; end
grid on;
set(gca,'XTick',0:20,'XTickLabel',d);
ylim([-50 400]);
line(xlim,[0 0]);
xlabel('Population Effect Size');
ylabel('Average Effect Size Inflation (% increase)');
legend('n = 20','n = 50','n = 100','n = 200','n = 1000','n = 2000');

%% Save plot
if ~isempty(fname2save)
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end
