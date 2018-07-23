%%  PRELIMINARIES
% -------------
% 0. Get a GitHub account fromhttps://github.com and send me your login name
% 1. Install GitHub client from https://desktop.github.com
% 2. Get repositories:
%       Psychtoolbox-3 - https://github.com/Psychtoolbox-3/Psychtoolbox-3
%       SilentSubstitutionToolbox - https://github.com/spitschan/SilentSubstitutionToolbox
%       MelanopicColourMatching - https://github.com/spitschan/MelanopicColourMatching
% 3. Add these repositories to the path

% READING
% -------
% http://cvrl.ioo.ucl.ac.uk/database/text/intros/introcmfs.htm
% http://color.psych.upenn.edu/brainard/papers/Brainard_Stockman_Colorimetry.pdf
%   10.3 Fundamentals of colorimetry
%   104. Color ordinate systems ? "Colour-Matching Functions" and "Cone Fundamentals"
% http://www.ling.upenn.edu/courses/ling525/color_vision.html


%% CONE SPECTRAL SENSITIVITIES
% ---------------------------

clear

% Wavelengths
wls = (390:1:830)'; 
S = WlsToS(wls);

% Change age
% receptorObj20 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 20);
% receptorObj30 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 30);
% receptorObj40 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 40);
% receptorObj50 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 50);
% receptorObj60 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 60);

% Plot one of them
% plot(wls, receptorObj20.T.T_energyNormalized', '--k')

% Change field size
% receptorObjFieldSize2Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 2);
% receptorObjFieldSize10Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);

% Plot one of them
% plot(wls, receptorObjFieldSize2Deg.T.T_energyNormalized', '--k');

% Export onto new variables for easy access
receptorObj = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10); % 10 degree field 
T_receptors = receptorObj.T.T_energyNormalized; % using normalized cone sensitivities
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

%% S&B COLOR MATCHING DATA
% -------------------
% If using Stiles&Burch color matching data, run this section of the
% code. If not, skip this section and use the next section of the code
% generates CMFs for a chosen set of primaries.

% Load S&B color matching data
mat = csvread('sbrgb10w.csv');
test_wls = mat(:, 1);
CMF = mat(:, 2:4);
r_match = round(645.16);
g_match = round(526.32);
b_match = round(444.44);
primaries = [r_match g_match b_match];


%% GENERATE CMF FOR A CHOSEN SET OF PRIMARIES
% -------------------------------------------
% If using Stiles&Burch color matching data, skip this section of the code.

% set up primaries
r_match = round(645); 
g_match = round(526); 
b_match = round(480);
primaries = [r_match g_match b_match];

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
    [CMF(i,:), T_fval(i)] = fmincon(@(CMF)opt_primaries_diff...
        (CMF,test_wls(i),primaries,[L;M;S],wls),startp,...
        [],[],[],[],minp,maxp,[], options);
end

% plot CMF
figure; hold on;
plot(test_wls,CMF(:,1), 'r', 'LineWidth', 2)
plot(test_wls,CMF(:,2), 'g', 'LineWidth', 2)
plot(test_wls,CMF(:,3), 'b', 'LineWidth', 2)

%% CONSTRUCT TEST AND MATCH SPECTRA
% Construct match and test spectra (assuming Gaussian shape)

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

% % if using delta function to construct spectra uncomment below:
% for i=1:n_test
%     neg_prim = find(CMF(i,:)<0);
%     pos_prim(:) = find(CMF(i,:)>0);
%     
%     spdt(wls2(:,i)==test_wls(i),i)=1;  % test light
%     spdt(wls2(:,i)==primaries(neg_prim),i)=CMF(i,neg_prim); % negative primary
%     
%     spdm(wls2(:,i)==primaries(pos_prim(1)),i)=CMF(i,pos_prim(1));  % positive primaries
%     spdm(wls2(:,i)==primaries(pos_prim(2)),i)=CMF(i,pos_prim(2));  
% end


% % plot test and match spectra
% figure; hold on; 
% subplot(2,1,1); plot(wls,spdt); title('test spectra'); pbaspect([3 1 1]);
% subplot(2,1,2); plot(wls,spdm); title('match spectra'); pbaspect([3 1 1]);
% 
% % plot test and match spectrum for a given wavelength
% wavelength = 500;  % choose wavelength
% wavelength_ind = find(test_wls==wavelength);
% figure; hold on; 
% subplot(2,1,1); hold on; plot(wls,spdt(:,wavelength_ind)); title('test spectrum'); 
% plot(xlim, [0,0], 'k:'); pbaspect([3 1 1]);
% subplot(2,1,2); hold on; plot(wls,spdm(:,wavelength_ind)); title('match spectrum'); 
% plot(xlim, [0,0], 'k:'); pbaspect([3 1 1]);


%% MELANOPSIN AND CONE STIMULATION
% --------------------------------

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


%% PLOT RECEPTOR STIMULATION AND CONTRAST
% ---------------------------------------

% plot cone and melanopic stimulation as a function of wavelength
figure; subplot(4,2,[1 2 3 4]); hold on; 
plot(test_wls,L_stimt,'r'); plot(test_wls,L_stimm,'rx');
plot(test_wls,M_stimt,'g'); plot(test_wls,M_stimm,'gx');  
plot(test_wls,S_stimt,'b'); plot(test_wls,S_stimm,'bx'); 
plot(test_wls(mel_wls),mel_stimt,'c'); plot(test_wls(mel_wls),mel_stimm,'cx'); 
legend({'L cone test', 'L cone macth', 'M cone test', 'M cone match', ...
    'S cone test', 'S cone match', 'melanopsin test', 'melanopsin match'})
title('receptor stimulation'); plot(xlim, [0,0], 'k:'); 
% zooming in on individuals receptor types
subplot(4,2,5); hold on; plot(test_wls,L_stimt,'r'); plot(test_wls,L_stimm,'rx'); 
title('L cone stimulation'); plot(xlim, [0,0], 'k:'); 
subplot(4,2,6); hold on; plot(test_wls,M_stimt,'g'); plot(test_wls,M_stimm,'gx'); 
title('M cone stimulation'); plot(xlim, [0,0], 'k:');
subplot(4,2,7); hold on; plot(test_wls,S_stimt,'b'); plot(test_wls,S_stimm,'bx'); 
title('S cone stimulation'); plot(xlim, [0,0], 'k:'); 
subplot(4,2,8); hold on; plot(test_wls(mel_wls),mel_stimt,'c'); 
plot(test_wls(mel_wls),mel_stimm,'cx'); 
title('melanopsin stimulation'); plot(xlim, [0,0], 'k:'); 

% % plot difference in receptor stimulation as a function of test wavelength
% figure; hold on; pbaspect([2 1 1]);
% plot(test_wls(mel_wls),mel_diff, 'c', 'LineWidth', 2)
% plot(test_wls,L_diff, 'r--')
% plot(test_wls,M_diff, 'g--')
% plot(test_wls,S_diff, 'b--')
% plot(xlim, [0,0], 'k:')
% legend({'melanopic diff', 'L cone diff', 'M cone diff', 'S cone diff', ''})

% % plot receptor contrast as a function of test wavelength
% figure; hold on; pbaspect([2 1 1]);
% plot(test_wls(mel_wls),mel_cont, 'c', 'LineWidth', 2)
% plot(test_wls,L_cont, 'r--')
% plot(test_wls,M_cont, 'g--')
% plot(test_wls,S_cont, 'b--')
% plot(xlim, [0,0], 'k:')
% legend({'melanopic contr', 'L cone contr', 'M cone contr', 'S cone contr', ''})

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

% % plot S cone contrast separately
% figure; hold on
% subplot(2,1,1); hold on; plot(test_wls,S_cont, 'b', 'LineWidth', 2); 
% title('S cone contrast'); plot(xlim, [0,0], 'k:')
% subplot(2,1,2); hold on; plot(test_wls,S_stimt, 'b'); 
% plot(test_wls,S_stimm, 'bx'); 
% plot(xlim, [0,0], 'k:'); title('S cone stimulation'); 
% legend({'S cone test', 'S cone match',''})




%% CMF GENERATION ALGORITHM FUNCTION
%-----------------------------------

% This section contains the function for the algorithm generating the CMF
% for a given set of non-overlapping primaries

function sqdiff = opt_primaries_diff(CMF,wls_test,primaries,cone_sens,range)

if sum(CMF<0)~=1
    sqdiff = inf;   % catch if more than one primary is negative
else
    
    % sort negative and positive primaries
    neg_prim = find(CMF<0);
    pos_prim(:) = find(CMF>0);
    
    wls=range;
    n_wls = length(wls);
    spdt = zeros(n_wls,1);  % test spectrum var
    spdm = zeros(n_wls,1);  % match spectrum var

    FWHM = 10; % half width of gaussian based on interference filter properties
    sigma = FWHM/2.4; 
    gauss = @(height,position)height.*exp(-((wls-position).^2)/(2*sigma^2));
      
    % construct test spectrum
    spdt = spdt+gauss(1,wls_test);
    spdt = spdt+gauss(abs(CMF(neg_prim)),primaries(neg_prim));

    % construct match spectrum
    spdm = spdm+gauss(CMF(pos_prim(1)),primaries(pos_prim(1)));
    spdm = spdm+gauss(CMF(pos_prim(2)),primaries(pos_prim(2)));
    
    % calculate cone stimulation
    stimm = cone_sens(1:3,:)*spdm;
    stimt = cone_sens(1:3,:)*spdt;

    % calculate sum of squared differences in stimulation
    sqdiff = sum((stimt-stimm).^2);
    
end
end
