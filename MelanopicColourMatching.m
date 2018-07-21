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
% wls = (380:1:780)';
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
receptorObjFieldSize10Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);

% Plot one of them
% plot(wls, receptorObjFieldSize2Deg.T.T_energyNormalized', '--k');

% Export onto new variables for easy access
T_receptors = receptorObjFieldSize10Deg.T.T_energy;
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

% 
% mat0 = csvread('linss10e_1.csv');
% mat0 = csvread('ss10e_1.csv');
 
% T_receptors = receptorObjFieldSize10Deg.T.T_energyNormalized;
% L = mat0(:,2)';
% M = mat0(:,3)';
% S = mat0(:,4)';
% mel = T_receptors(4,:);


% COLOR MATCHING DATA
% -------------------
mat = csvread('sbrgb10w.csv');
wls_snb = mat(:, 1);
snb_rgb = mat(:, 2:4);


% Construct match and test spectra
r_snb = round(645.16);
g_snb = round(526.32);
b_snb = round(444.44);
rgb_wls = [r_snb g_snb b_snb];

n_snb = length(wls_snb);
n_wls = length(wls);
wls2 = repmat(wls,1,n_snb);
spdt = zeros(n_wls,n_snb);
spdm = zeros(n_wls,n_snb);

% constructing spectra assuming a gaussian shape of primaries and test
FWHM = 10; % FWMH of gaussian based on interference filter properties
sigma = FWHM/2.4; % stdev of gaussian
gauss = @(height,position)height.*exp(-((wls-position).^2)/(2*sigma^2));

for i=1:n_snb
    neg_prim = find(snb_rgb(i,:)<0);
    pos_prim(:) = find(snb_rgb(i,:)>0);
    
    spdt(:,i) = spdt(:,i)+gauss(1,wls_snb(i));
    spdt(:,i) = spdt(:,i)+gauss(abs(snb_rgb(i,neg_prim)),rgb_wls(neg_prim));
    
    spdm(:,i) = spdm(:,i)+gauss(snb_rgb(i,pos_prim(1)),rgb_wls(pos_prim(1)));
    spdm(:,i) = spdm(:,i)+gauss(snb_rgb(i,pos_prim(2)),rgb_wls(pos_prim(2)));
end


% % if using delta function to construct spectra uncomment below:
% for i=1:n_snb
%     spd1(wls2(:,i)==wls_snb(i),i)=1;  % test light
%     spd2(wls2(:,i)==r_snb,i)=snb_rgb(i,1);  % red priamry
%     spd2(wls2(:,i)==g_snb,i)=snb_rgb(i,2);  % green primary
%     spd2(wls2(:,i)==b_snb,i)=snb_rgb(i,3);  % blue primary
% end

% plot test and match spectra
figure; hold on; 
subplot(2,1,1); plot(wls,spdt); title('test spectra'); pbaspect([3 1 1]);
subplot(2,1,2); plot(wls,spdm); title('match spectra'); pbaspect([3 1 1]);
% plot one example test and match spectrum
figure; hold on; 
subplot(2,1,1); hold on; plot(wls,spdt(:,52)); title('test spectrum'); 
plot(xlim, [0,0], 'k:'); pbaspect([3 1 1]);
subplot(2,1,2); hold on; plot(wls,spdm(:,52)); title('match spectrum'); 
plot(xlim, [0,0], 'k:'); pbaspect([3 1 1]);


% MELANOPSIN AND CONE STIMULATION
% -------------------

% calculate cone and melanopsin stimulation
L_stimt = L*spdt;
L_stimm = L*spdm;
M_stimt = M*spdt;
M_stimm = M*spdm;
S_stimt = S*spdt;
S_stimm = S*spdm;
mel_stimt = mel*spdt;
mel_stimm = mel*spdm;
% plot cone and melanopic stimulation as a function of wavelength
figure; hold on; 
plot(wls_snb,L_stimt,'r'); plot(wls_snb,L_stimm,'rx');
plot(wls_snb,M_stimt,'g'); plot(wls_snb,M_stimm,'gx');  
plot(wls_snb,S_stimt,'b'); plot(wls_snb,S_stimm,'bx'); 
plot(wls_snb,mel_stimt,'c'); plot(wls_snb,mel_stimm,'cx'); 
legend({'L cone test', 'L cone macth', 'M cone test', 'M cone match', ...
    'S cone test', 'S cone match', 'melanopsin test', 'melanopsin match'})
% zooming in on individuals receptor types
figure; hold on; 
subplot(2,2,1); hold on; plot(wls_snb,L_stimt,'r'); plot(wls_snb,L_stimm,'rx'); 
title('L cone stimulation'); plot(xlim, [0,0], 'k:'); pbaspect([2 1 1]);
subplot(2,2,2); hold on; plot(wls_snb,M_stimt,'g'); plot(wls_snb,M_stimm,'gx'); 
title('M cone stimulation'); plot(xlim, [0,0], 'k:'); pbaspect([2 1 1]);
subplot(2,2,3); hold on; plot(wls_snb,S_stimt,'b'); plot(wls_snb,S_stimm,'bx'); 
title('S cone stimulation'); plot(xlim, [0,0], 'k:'); pbaspect([2 1 1]);
subplot(2,2,4); hold on; plot(wls_snb,mel_stimt,'c'); plot(wls_snb,mel_stimm,'cx'); 
title('melanopsin stimulation'); plot(xlim, [0,0], 'k:'); pbaspect([2 1 1]);

% calculate cone and melanopsin contrast
L_cont = (L_stimt-L_stimm)./(L_stimm);
M_cont = (M_stimt-M_stimm)./(M_stimm);
S_cont = (S_stimt-S_stimm)./(S_stimm);
mel_cont = (mel_stimt-mel_stimm)./(mel_stimm);

% plot contrast as a function of test wavelength
figure; hold on; pbaspect([2 1 1]);
plot(wls_snb,mel_cont, 'c', 'LineWidth', 2)
plot(wls_snb,L_cont, 'r--')
plot(wls_snb,M_cont, 'g--')
plot(wls_snb,S_cont, 'b--')
plot(xlim, [0,0], 'k:')
legend({'melanopic contr', 'L cone contr', 'M cone contr', 'S cone contr', ''})

% plot S cone and melanopsin contrast separately
figure; hold on
subplot(2,1,1); hold on; plot(wls_snb,mel_cont, 'c', 'LineWidth', 2); 
title('Melanopic contrast'); plot(xlim, [0,0], 'k:')
subplot(2,1,2); hold on; plot(wls_snb,mel_stimt, 'c'); 
plot(wls_snb,mel_stimm, 'cx'); 
plot(xlim, [0,0], 'k:');  title('Melanopsin stimulation'); 
legend({'melanopsin test', 'melanopsin match',''})

figure; hold on
subplot(2,1,1); hold on; plot(wls_snb,S_cont, 'b', 'LineWidth', 2); 
title('S cone contrast'); plot(xlim, [0,0], 'k:')
subplot(2,1,2); hold on; plot(wls_snb,S_stimt, 'b'); 
plot(wls_snb,S_stimm, 'bx'); 
plot(xlim, [0,0], 'k:'); title('S cone stimulation'); 
legend({'S cone test', 'S cone match',''})

% calculate cone and melanopsin difference in stimulation
L_diff = (L_stimt-L_stimm);
M_diff = (M_stimt-M_stimm);
S_diff = (S_stimt-S_stimm);
mel_diff = (mel_stimt-mel_stimm);

% plot difference as a function of test wavelength
figure; hold on; pbaspect([2 1 1]);
plot(wls_snb,mel_diff, 'c', 'LineWidth', 2)
plot(wls_snb,L_diff, 'r--')
plot(wls_snb,M_diff, 'g--')
plot(wls_snb,S_diff, 'b--')
plot(xlim, [0,0], 'k:')
legend({'melanopic diff', 'L cone diff', 'M cone diff', 'S cone diff', ''})




%%

% CHNAGING PRIMARIES
% -------------------

clear 

wls = (390:1:830)';
S = WlsToS(wls);
receptorObjFieldSize10Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);
T_receptors = receptorObjFieldSize10Deg.T.T_energyNormalized;
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

r_match = round(645.16);
g_match = round(526.32);
b_match = round(444.44);
primaries = [r_match,g_match,b_match];

wls_test = (390:5:830)';
n_test = length(wls_test);
    
startp = [-1,1,1; 1,-1,1; 1,1,-1];
minp = [-5,-5,-5];
maxp = [5,5,5];
% options = optimset('TolCon', 1e-10000);
options = optimoptions('fmincon','Display',        'none',...
    'OptimalityTolerance',      0,...
    'StepTolerance',            0,...
    'MaxIterations',            3e3,...
    'MaxFunctionEvaluations',   3e3);

tic                   

% options = optimoptions('fmincon', 'Display','iter','Diagnostics','on');

for i=1:n_test
    for c = 1:3
        [return_rgb(i,:,c),fval(i,c)] = fmincon(@(T_rgb)opt_primaries...
            (T_rgb,wls_test(i),primaries,[L;M;S]),startp(c,:),...
            [],[],[],[],minp,maxp,[], options);
    end
    [sumsqcontr(i),position(i)] = min(fval(i,:));
    T_rgb(i,:) = return_rgb(i,:,position(i));
end



% for l=1:n_test
%     [L_cont(l) M_cont(l) S_cont(l)] = opt_primaries_return(T_rgb(l,:),wls_test(l),primaries,[L;M;S]);
% end

toc
%% 

figure
r = plot(wls_test,snb_rgb(:,1),'r','LineWidth',5); hold on
rm = plot(wls_test,T_rgb(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(wls_test,snb_rgb(:,2),'g','LineWidth',5); hold on
gm = plot(wls_test,T_rgb(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(wls_test,snb_rgb(:,3),'b','LineWidth',5); hold on
bm = plot(wls_test,T_rgb(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;
title('using sum of squared (non-normalized) differences')
legend({'r SNB','r model','g SNB','g model','b SNB','b model'})

%%

function sqdiff = opt_primaries(T_rgb,wls_test,primaries,cone_sens)

if sum(T_rgb<0)~=1
    sqdiff = inf;    
else
    wls = (390:1:830)';

    n_wls = length(wls);
    spdt = zeros(n_wls,1);
    spdm = zeros(n_wls,1);

    FWHM = 10; 
    sigma = FWHM/2.4; 
    gauss = @(height,position)height.*exp(-((wls-position).^2)/(2*sigma^2));
        
    neg_prim = find(T_rgb<0);
    pos_prim(:) = find(T_rgb>0);

    spdt = spdt+gauss(1,wls_test);
    spdt = spdt+gauss(abs(T_rgb(neg_prim)),primaries(neg_prim));

    spdm = spdm+gauss(T_rgb(pos_prim(1)),primaries(pos_prim(1)));
    spdm = spdm+gauss(T_rgb(pos_prim(2)),primaries(pos_prim(2)));


    stimm = cone_sens(1:3,:)*spdm;
    stimt = cone_sens(1:3,:)*spdt;

%     sqdiff = sum((stimt-stimm).^2);
    
    L = cone_sens(1,:);
    M = cone_sens(2,:);
    S = cone_sens(3,:);
    L_stimt = L*spdt;
    L_stimm = L*spdm;
    M_stimt = M*spdt;
    M_stimm = M*spdm;
    S_stimt = S*spdt;
    S_stimm = S*spdm;
    L_cont = ((L_stimt-L_stimm)./(L_stimm))*100;
    M_cont = ((M_stimt-M_stimm)./(M_stimm))*100;
    S_cont = ((S_stimt-S_stimm)./(S_stimm))*100;
    
    sqdiff = (L_cont.^2+M_cont.^2+S_cont.^2);
end
end

