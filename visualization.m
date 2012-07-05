% this script simply takes the saved data and creates the exact plots as found in the paper and thesis
% Remco Tukker

load data2.mat

%%
%ch = ranksum(squeeze(perfp(1,6,1:50)),squeeze(perfp(2,6,1:50)), 'method', 'exact')
chp3 = ranksum(squeeze(perfp(1,3,:)),squeeze(perfp(2,3,:)))  %note, this is a normal approximation; exact method takes too long for 200 samples (50 is doable)
chp6 = ranksum(squeeze(perfp(1,6,:)),squeeze(perfp(2,6,:)))

cht4 = ranksum(squeeze(perft(1,4,:)),squeeze(perft(2,4,:)))
cht10 = ranksum(squeeze(perft(1,10,:)),squeeze(perft(2,10,:)))

%% temp visualization


 cmat = [0.3 0.3 1; 1 0.5 0];
 grouping2 = {'      0' '      0' '    1' '    1' '    2' '    2' '    3' '    3' '    4' '    4' '    5' ...
     '    5' '    6' '    6' '    7' '    7' '    8' '    8' '    15' '    15' '    30' '    30' '    90' '    90'};
 g3 = {'u','y'}; 
 grouping3 = [g3 g3 g3 g3 g3 g3 g3 g3 g3 g3 g3 g3];
 
 figure('Position',[1 1 780 600],'Color','w')
 boxplot(perfcollectiont', { grouping2  grouping3},'plotstyle' ,'compact', 'notch', 'off', 'factorgap' , [5 0], 'colorgroup', grouping3, ...
     'labels', grouping2, 'colors', cmat , 'symbol', '.', 'labelverbosity', 'minor', 'outliersize', 2, 'medianstyle','line', ...
     'jitter', 0, 'labelorientation', 'horizontal')
xlabel('Memory Length','fontsize',13);
ylabel('Performance','fontsize',13);
set(gca,'YGrid','on','Ycolor',[0.4 0.4 0.4]); %grey grid
ylim([0.35 1]);
box off;
set(findobj(gca,'Type','text'),'FontSize',12)
set(gca,'FontSize',12);

 % make median lines black and big
set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',3);
set(findobj(gcf,'Tag','Box'),'LineWidth',8);
set(findobj(gcf,'Tag','Whisker'),'LineWidth',2);

bars = findobj(gcf,'Tag','Box');
legend([bars(1), bars(2)], 'Uniform Input' , 'Topological Input', 'FontSize',13);

 % make outlier dots gray and big
set(findobj(gcf,'Tag','Outliers'),'MarkerSize',5);

Caxes = copyobj(gca,gcf); %make the rest of the stuff black again
xlabel(Caxes,'');
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 

line([0 100],[.5 .5],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line


%bar
figure('Position',[1 1 780 600],'Color','w')

bar(barperft, 'FaceColor', [0.3 0.3 1])
set(gca,'xticklabel', {0 1 2 3 4 5 6 7 8 15 30 90} );
xlabel('Memory Length','fontsize',13);
ylabel('Difference in Median Performance (Topological - Uniform)','fontsize',13);
set(gca,'FontSize',12);

box off;
%set(gcf, );
%ylim([0 1]);
set(gca,'YGrid','on','Ycolor',[0.4 0.4 0.4]); %grey grid

Caxes = copyobj(gca,gcf); %make the rest of the stuff black again
xlabel(Caxes,'');
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 


%% pol visualization


 cmat = [0.3 0.3 1; 1 0.5 0];
 grouping2 = {'          1' '          1' '        2' '        2' '         3' '         3' '         4' '         4' '        5' '        5' ...
     '      6' '      6' '     7' '     7' '     8' '     8' }; 
     %'    5' '    5' '    6' '    6' '    7' '    7' '    8' '    8' '    15' '    15' '    30' '    30'};
 g3 = {'u','y'}; 
 grouping3 = [g3 g3 g3 g3 g3 g3 g3 g3];
 
 figure('Position',[1 1 780 600],'Color','w')
 boxplot(perfcollectionp', { grouping2  grouping3},'plotstyle' ,'compact', 'notch', 'off', 'factorgap' , [5 0], 'colorgroup', grouping3, ...
     'labels', grouping2, 'colors', cmat , 'symbol', '.', 'labelverbosity', 'minor', 'outliersize', 2, 'medianstyle','line', ...
     'jitter', 0, 'labelorientation', 'horizontal')
xlabel('Number of Low-Level Cues','fontsize',13);
ylabel('Performance','fontsize',13);
set(gca,'YGrid','on','Ycolor',[0.4 0.4 0.4]); %grey grid
ylim([0.0 1]);
box off;
set(findobj(gca,'Type','text'),'FontSize',12)
set(gca,'FontSize',12);

 % make median lines black and big
set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',3);
set(findobj(gcf,'Tag','Box'),'LineWidth',12);
set(findobj(gcf,'Tag','Whisker'),'LineWidth',2);

bars = findobj(gcf,'Tag','Box');  %whisker
legend([bars(1), bars(2)], 'Uniform Input' , 'Topological Input', 'FontSize',13);

 % make outlier dots gray and big
set(findobj(gcf,'Tag','Outliers'),'MarkerSize',5);

Caxes = copyobj(gca,gcf); %make the rest of the stuff black again
xlabel(Caxes,'');
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 

line([0 100],[.5 .5],[-1 -1],'color',[255/255, 95/255, 121/255],'linewidth',1); %chance level line

line([0.55 0.6+21.15/8],[.25 .25],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8 0.6+2*21.15/8],[.125 .125],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8*2 0.6+21.15/8*3],[.0625 .0625],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8*3 0.6+21.15/8*4],[1/2^5 1/2^5],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8*4 0.6+21.15/8*5],[1/2^6 1/2^6],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8*5 0.6+21.15/8*6],[1/2^7 1/2^7],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8*6 0.6+21.15/8*7],[1/2^8 1/2^8],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line
line([0.6+21.15/8*7 0.6+21.15/8*8],[1/2^9 1/2^9],[-1 -1],'color',[215/255, 25/255, 28/255],'linewidth',2); %chance level line


%bar
figure('Position',[1 1 780 600],'Color','w')

bar(barperfp, 'FaceColor', [0.3 0.3 1])
set(gca,'xticklabel',{1 2 3 4 5 6 7 8});
xlabel('Number of Low-Level Cues','fontsize',13);
ylabel('Difference in Median Performance (Topological - Uniform)','fontsize',13);
set(gca,'FontSize',12);

box off;
%set(gcf, );
%ylim([0 1]);
set(gca,'YGrid','on','Ycolor',[0.4 0.4 0.4]); %grey grid

Caxes = copyobj(gca,gcf); %make the rest of the stuff black again
xlabel(Caxes,'');
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 

