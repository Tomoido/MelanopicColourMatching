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


% CONE SPECTRAL SENSITIVITIES
% ---------------------------
wls = (380:1:780)';
S = WlsToS(wls);

% Change age
receptorObj20 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 20);
receptorObj30 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 30);
receptorObj40 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 40);
receptorObj50 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 50);
receptorObj60 = SSTReceptorHuman('S', S, 'obsAgeInYrs', 60);

% Plot one of them
plot(wls, receptorObj20.T.T_energyNormalized', '--k')

% Change field size
receptorObjFieldSize2Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 2);
receptorObjFieldSize10Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);

% Plot one of them
plot(wls, receptorObjFieldSize2Deg.T.T_energyNormalized', '--k');