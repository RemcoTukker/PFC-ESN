% This file calls ESNtest4 a large number of times with a couple of settings. This gives the performance distributions found in the thesis and paper
% Created by Remco Tukker

clear all;

trials = 200;   %number of trials per parameter setting

%alternative parameter settings
%p = [0.8 200 0.05 0.1 0.9 5 0.3 15 0.25]; %better for p
%p = [1.0 200 0.05 0.1 0.9 5 0.3 15 0.25]; %better for t


%% policy abstraction

p = [0.6 200 0.05 0.1 0.9 5 0.3 15 0.25];

nrofcues =   [1 2 3 4 5 6 7 8];
nrofinputs = [3 5 7 9 11 13 15 17];

% the performance, mean output weight, and difference between training and test performance is returned
perfp = zeros(2, numel(nrofcues), trials);
meansp = zeros(2, numel(nrofcues), trials);
diffsp = zeros(2, numel(nrofcues), trials);


for i = 1:trials
   for cues = 1:numel(nrofcues)
       [perfp(1, cues, i) meansp(1, cues, i) diffsp(1, cues, i) ] =  ESNtest4(1, 0, 1, nrofcues(cues), nrofinputs(cues), p(1), p(2), p(3), p(4) , p(5), p(6), p(7), p(8), p(9));
       [perfp(2, cues, i) meansp(2, cues, i) diffsp(2, cues, i) ] =  ESNtest4(1, 0, 2, nrofcues(cues), nrofinputs(cues), p(1), p(2), p(3), p(4) , p(5), p(6), p(7), p(8), p(9));
   end
   disp(i);
end

%put the data in some nice collection matrices: perfp dimension 1 is the difference in input method, dim 2 is the number of cues and dim 3 are the individual trials
perfcollectionp = zeros(2* size(perfp,2), trials );
for i = 1:size(perfp,2)
   perfcollectionp(i*2-1,:) = perfp(1,i,:);
   perfcollectionp(i*2,:) = perfp(2,i,:);
end

% take the median as a measure of performance difference
barperfp = median(squeeze(perfp(1,:,:)),2) - median(squeeze(perfp(2,:,:)),2);

% %make picture with current parameter settings somewhere
% figure(j)
% a = subplot(2,1,1);
% boxplot(perfcollectionp', 'notch', 'on', 'colorgroup', repmat([0 1], 1, numel(nrofcues)) )
% b = subplot(2,1,2);
% bar(barperfp)
% title(a, ['TDs ',num2str(p(j,1)),' n ',num2str(p(j,2)) , ' c ',num2str(p(j,3)),' ic ',num2str(p(j,4)), ' r ',num2str(p(j,5)), ...
%         ' is ',num2str(p(j,6)), ' bs ',num2str(p(j,7)), ' ts ',num2str(p(j,8)) , ' tr ' ,num2str(p(j,9))  ])    
    
%% temporal abstraction

% (exactly the same structure as above; see the comments there)

p = [1.0 100 0.05 0.1 0.9 5 0.3 15 0.25];

memlengths = [0 1 2 3 4 5 6 7 8 15 30 90];

perft = zeros(2, numel(memlengths), trials);
meanst = zeros(2, numel(memlengths), trials);
diffst = zeros(2, numel(memlengths), trials);

for i = 1:trials 
    
   for meml = 1:numel(memlengths)
       [perft(1, meml, i) meanst(1, meml, i) diffst(1, meml, i) ] = ESNtest4(2, memlengths(meml), 1, 3, 3, p(1), p(2), p(3), p(4) , p(5), p(6), p(7), p(8), p(9) );
       [perft(2, meml, i) meanst(2, meml, i) diffst(2, meml, i) ] = ESNtest4(2, memlengths(meml), 2, 3, 3, p(1), p(2), p(3), p(4) , p(5), p(6), p(7), p(8), p(9) );
   end
   
   disp(i)

end

perfcollectiont = zeros(2* size(perft,2), trials );
for i = 1:size(perft,2)
    perfcollectiont(i*2-1,:) = perft(1,i,:);
    perfcollectiont(i*2,:) = perft(2,i,:);
end

barperft = median(squeeze(perft(1,:,:)),2) - median(squeeze(perft(2,:,:)),2);
% 
% figure(1)
% boxplot(perfcollectiont', 'notch', 'on', 'colorgroup', repmat([0 1], 1, numel(memlengths)) )
% 
% figure(2)
% bar(barperft)

%% 

save data2;
