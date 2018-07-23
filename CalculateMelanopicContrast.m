% MELANOPIC CONTRAST ACROSS CMFS
%-------------------------------

% This function returns the melanopic contrast for the CMF of a given set
% of priamries. It takes a three element vector of primaries (rgb in nm) as input. 

function mel_cont = CalculateMelanopicContrast(primaries)

% Get cone sensitivities
wls = (390:1:830)'; 
S = WlsToS(wls);
receptorObj = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10); % 10 degree field 
T_receptors = receptorObj.T.T_energyNormalized; % using normalized cone sensitivities
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

% Generate CMF for chosen primaries
test_wls = (390:5:810)';
n_test = length(test_wls);
startp = [-1,1,1];
minp = [-10,-10,-10];
maxp = [20,20,20];

options = optimoptions('fmincon','Display', 'iter',...
    'Algorithm','sqp', ...
    'OptimalityTolerance',      0,...
    'StepTolerance',            0,...
    'MaxIterations',            3e3,...
    'MaxFunctionEvaluations',   3e3);

for i=1:n_test
    [CMF(i,:), T_fval(i)] = fmincon(@(CMF)GenerateCMF...
        (CMF,test_wls(i),primaries,[L;M;S],wls),startp,...
        [],[],[],[],minp,maxp,[], options);
end

% plot CMF
figure; hold on;
plot(test_wls,CMF(:,1), 'r', 'LineWidth', 2)
plot(test_wls,CMF(:,2), 'g', 'LineWidth', 2)
plot(test_wls,CMF(:,3), 'b', 'LineWidth', 2)

% construct test and match spectra
n_test = length(test_wls);
n_wls = length(wls);
wls2 = repmat(wls,1,n_test);
spdt = zeros(n_wls,n_test);
spdm = zeros(n_wls,n_test);

FWHM = 10; % FWMH of gaussian based on interference filter properties
sigma = FWHM/2.4; % stdev of gaussian
gauss = @(height,position)height.*exp(-((wls-position).^2)/(2*sigma^2));

for i=1:n_test
    neg_prim = find(CMF(i,:)<0);  % identify negative primary
    pos_prim(:) = find(CMF(i,:)>0);  % identify positive primaries
    
    % test spectrum using test light and negative primary
    spdt(:,i) = spdt(:,i)+gauss(1,test_wls(i));
    spdt(:,i) = spdt(:,i)+gauss(abs(CMF(i,neg_prim)),primaries(neg_prim));
    
    % match spectrum using positive primaries
    spdm(:,i) = spdm(:,i)+gauss(CMF(i,pos_prim(1)),primaries(pos_prim(1)));
    spdm(:,i) = spdm(:,i)+gauss(CMF(i,pos_prim(2)),primaries(pos_prim(2)));
end

% calculate cone and melanopsin stimulation
L_stimt = L*spdt;
L_stimm = L*spdm;
M_stimt = M*spdt;
M_stimm = M*spdm;
S_stimt = S*spdt;
S_stimm = S*spdm;
mel_stimt = mel*spdt;
mel_stimm = mel*spdm;
mel_wls = find(test_wls==390):find(test_wls==620);  % not calculating melanopsin 
% stimulation 3 log units above and below the nominal max sensitivity
mel_stimt = mel_stimt(mel_wls);                  
mel_stimm = mel_stimm(mel_wls);

% calculate cone and melanopsin difference in stimulation
L_diff = (L_stimt-L_stimm);
M_diff = (M_stimt-M_stimm);
S_diff = (S_stimt-S_stimm);
mel_diff = (mel_stimt-mel_stimm);

% calculate cone and melanopsin contrast
L_cont = (L_stimt-L_stimm)./(L_stimm);
M_cont = (M_stimt-M_stimm)./(M_stimm);
S_cont = (S_stimt-S_stimm)./(S_stimm);
mel_cont = (mel_stimt-mel_stimm)./(mel_stimm);

% plot cone and melanopic stimulation as a function of wavelength
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = ['For primaries: r(' num2str(round(primaries(1))) ') g('...
    num2str(round(primaries(2))) ') b(' num2str(round(primaries(3))) ')']; 
p.TitlePosition = 'centertop';  p.FontSize = 14; p.FontWeight = 'bold';
hold on;  
subplot(4,2,[1 2 3 4],'Parent',p); hold on; 
plot(test_wls,L_stimt,'r'); plot(test_wls,L_stimm,'rx');
plot(test_wls,M_stimt,'g'); plot(test_wls,M_stimm,'gx');  
plot(test_wls,S_stimt,'b'); plot(test_wls,S_stimm,'bx'); 
plot(test_wls(mel_wls),mel_stimt,'c'); plot(test_wls(mel_wls),mel_stimm,'cx'); 
legend({'L cone test', 'L cone macth', 'M cone test', 'M cone match', ...
    'S cone test', 'S cone match', 'melanopsin test', 'melanopsin match'})
title('receptor stimulation'); plot(xlim, [0,0], 'k:'); 
% zooming in on individuals receptor types
subplot(4,2,5,'Parent',p); hold on; 
plot(test_wls,L_stimt,'r'); plot(test_wls,L_stimm,'rx'); 
title('L cone stimulation'); plot(xlim, [0,0], 'k:'); 
subplot(4,2,6,'Parent',p); 
hold on; plot(test_wls,M_stimt,'g'); plot(test_wls,M_stimm,'gx'); 
title('M cone stimulation'); plot(xlim, [0,0], 'k:');
subplot(4,2,7,'Parent',p); hold on; 
plot(test_wls,S_stimt,'b'); plot(test_wls,S_stimm,'bx'); 
title('S cone stimulation'); plot(xlim, [0,0], 'k:'); 
subplot(4,2,8,'Parent',p); hold on; 
plot(test_wls(mel_wls),mel_stimt,'c'); plot(test_wls(mel_wls),mel_stimm,'cx'); 
title('melanopsin stimulation'); plot(xlim, [0,0], 'k:'); 

% plot melanopsin contrast and difference in stimulation separately
figure; hold on
subplot(3,1,1); hold on; plot(test_wls(mel_wls),mel_stimt, 'c', 'LineWidth', 2); 
plot(test_wls(mel_wls),mel_stimm, 'cx', 'LineWidth', 2); 
plot(xlim, [0,0], 'k:');  title('Melanopsin stimulation'); 
legend({'melanopsin test', 'melanopsin match',''})
subplot(3,1,2); hold on; plot(test_wls(mel_wls),mel_diff, 'c', 'LineWidth', 2); 
title('Melanopic stimulation difference'); plot(xlim, [0,0], 'k:')
subplot(3,1,3); hold on; plot(test_wls(mel_wls),mel_cont, 'c', 'LineWidth', 2); 
title('Melanopic contrast'); plot(xlim, [0,0], 'k:')

