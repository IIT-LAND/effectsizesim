function effect_size_inflation_sim(fname2save)
%
%   fname2save = filename for pdf plot to save (leave empty if you don't want to save)
%

d = 0:0.1:2;
for i = 1:length(d)
    Results = effect_size_simulation(d(i));
    avg_es_inf(:,i) = Results.avg_es_inflate_percentincrease;
    close all; % close figure created by effect_size_simulation
end

lineWidth = 2;
% fontSize = 10;
% fontWeight = 'b';

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

if ~isempty(fname2save)
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end
