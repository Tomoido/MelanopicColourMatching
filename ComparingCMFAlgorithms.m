%% Comparison between 3 algorithms for CMF generation

% This code generates and compares CMF generation algorithms against Stiles
% and Burch experimental CMF data.

% This contains three algorithms for generation of CMFs:
% 1. Minimizing squared difference in normalized cone stimulation
% 2. Minimizing squared difference in absolute cone stimulation
% 3. Minimizing squared cone contrast


%% Run 3 CMF algorithms

clear 
mat = csvread('sbrgb10w.csv'); % load Stiles and Burch benchmark data
test_wls = mat(:, 1);
snb_CMF = mat(:, 2:4);
r_match = round(645.16);
g_match = round(526.32);
b_match = round(444.44);
primaries = [r_match,g_match,b_match];

% load cone sensitivities
wls = (390:1:830)'; 
S = WlsToS(wls);
receptorObj= SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);
T_receptors = receptorObj.T.T_energyNormalized; % starting with normlized
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);

% set up parameters for minimization algorithms
n_test = length(test_wls);
startp = [-1,1,1]; 
minp = [-5,-5,-5]; 
maxp = [5,5,5];
options = optimoptions('fmincon','Display', 'none',...
    'Algorithm','sqp', ...
    'OptimalityTolerance',      0,...
    'StepTolerance',            0,...
    'MaxIterations',            3e3,...
    'MaxFunctionEvaluations',   3e3);


% 1. Algorithm minimizing squared difference in normalized stimulation
for i=1:n_test
    [ndiff_CMF(i,:), ndiff_fval(i)] = fmincon(@(T_rgb)opt_primaries_diff...
        (T_rgb,test_wls(i),primaries,[L;M;S]),startp,...
        [],[],[],[],minp,maxp,[], options);
end


% 2. Algorithm minimizing squared difference in absolute stimulation
% loading non-normalized cone sensitivities
T_receptors = receptorObj.T.T_energy;
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

for i=1:n_test
    [diff_CMF(i,:), diff_fval(i)] = fmincon(@(T_rgb)opt_primaries_diff...
        (T_rgb,test_wls(i),primaries,[L;M;S]),startp,...
        [],[],[],[],minp,maxp,[], options);
end


% 3. Algorithm minimizing squared cone contrast
% adding more starting values: where each primary is negative once
startp = [-1,1,1; 1,-1,1; 1,1,-1];
for i=1:n_test
    for c = 1:3 % we go through this extra loop to load the extra starting values
        % this ensures that the algorithm won't get stuck in local minima
        % and will converge in a reasonable number of
        % iterations/f evaluations
        [return_CMF(i,:,c),fval(i,c)] = fmincon(@(T_rgb)opt_primaries_cont...
            (T_rgb,test_wls(i),primaries,[L;M;S]),startp(c,:),...
            [],[],[],[],minp,maxp,[], options);
    end
    [cont_fval(i),position(i)] = min(fval(i,:));
    cont_CMF(i,:) = return_CMF(i,:,position(i));
end

% saving output CMFs
save('algorithm_CMFs','cont_CMF','diff_CMF','ndiff_CMF')


%% Plot output CMFs against S&B data

figure

subplot(3,1,1)
r = plot(test_wls,snb_CMF(:,1),'r','LineWidth',5); hold on
rm = plot(test_wls,cont_CMF(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(test_wls,snb_CMF(:,2),'g','LineWidth',5); hold on
gm = plot(test_wls,cont_CMF(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(test_wls,snb_CMF(:,3),'b','LineWidth',5); hold on
bm = plot(test_wls,cont_CMF(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;
legend({'r SNB','r model','g SNB','g model','b SNB','b model'})
title('using sum of squared cone contrast')

subplot(3,1,2)
r = plot(test_wls,snb_CMF(:,1),'r','LineWidth',5); hold on
rm = plot(test_wls,diff_CMF(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(test_wls,snb_CMF(:,2),'g','LineWidth',5); hold on
gm = plot(test_wls,diff_CMF(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(test_wls,snb_CMF(:,3),'b','LineWidth',5); hold on
bm = plot(test_wls,diff_CMF(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;
title('using sum of squared differences (non-normalized)')

subplot(3,1,3)
r = plot(test_wls,snb_CMF(:,1),'r','LineWidth',5); hold on
rm = plot(test_wls,ndiff_CMF(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(test_wls,snb_CMF(:,2),'g','LineWidth',5); hold on
gm = plot(test_wls,ndiff_CMF(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(test_wls,snb_CMF(:,3),'b','LineWidth',5); hold on
bm = plot(test_wls,ndiff_CMF(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;
title('using sum of squared differences (normalized)')

%% Function for the difference minimization algorithms

function sqdiff = opt_primaries_diff(CMF,wls_test,primaries,cone_sens)

if sum(CMF<0)~=1
    sqdiff = inf;   % catch if more than one primary is negative
else
    
    % sort negative and positive primaries
    neg_prim = find(CMF<0);
    pos_prim(:) = find(CMF>0);
    
    wls = (390:1:830)';  % wavelengths for spectra
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

%% Function for the contrast minimization algorithms

function sqcont = opt_primaries_cont(CMF,wls_test,primaries,cone_sens)

if sum(CMF<0)~=1
    sqcont = inf;   % catch if more than one primary is negative    
else
    
    % sort negative and positive primaries
    neg_prim = find(CMF<0);
    pos_prim(:) = find(CMF>0);
    
    wls = (390:1:830)';  % wavelengths for spectra
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
    
    % calculate cone stimulation for test and match lights
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
    
    % calculate contrast
    sqcont = (L_cont.^2+M_cont.^2+S_cont.^2);
end
end
