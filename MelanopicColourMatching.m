%
% PRELIMINARIES
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


% CONE SPECTRAL SENSITIVITIES
% ---------------------------
wls = (380:1:780)';
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
receptorObjFieldSize10Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);

% Plot one of them
% plot(wls, receptorObjFieldSize2Deg.T.T_energyNormalized', '--k');

% Export onto new variables for easy access
T_receptors = receptorObjFieldSize10Deg.T.T_energy;
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

% COLOR MATCHING DATA
% -------------------
mat = csvread('sbrgb10w.csv');
wls_snb = mat(:, 1);
T_rgb = mat(:, 2:4);

% Construct match and test spectra
r_snb = round(645.16);
g_snb = round(526.32);
b_snb = round(444.44);

n_snb = length(wls_snb);
n_wls = length(wls);
wls2 = repmat(wls,1,n_snb);
spd1 = zeros(n_wls,n_snb);
spd2 = zeros(n_wls,n_snb);

for i=1:n_snb
    spd1(wls2(:,i)==wls_snb(i),i)=1;  % should this be 1?
    spd2(wls2(:,i)==r_snb,i)=T_rgb(i,1);
    spd2(wls2(:,i)==g_snb,i)=T_rgb(i,2);
    spd2(wls2(:,i)==b_snb,i)=T_rgb(i,3);
end

% Look at cone and melanopsin contrast
mel_cont = mel*(spd1-spd2)./(mel*(spd2));
L_cont = L*(spd1-spd2)./(L*(spd2));
M_cont = M*(spd1-spd2)./(M*(spd2));
S_cont = S*(spd1-spd2)./(S*(spd2));

% Plot as a function of contrast
figure; hold on;
plot(wls_snb,mel_cont, 'c', 'LineWidth', 2)
plot(wls_snb,L_cont, 'r--')
plot(wls_snb,M_cont, 'g--')
plot(wls_snb,S_cont, 'b--')
legend({'melanopic contr', 'L cone contr', 'M cone contr', 'S cone contr'})
plot(xlim, [0,0], 'k:')