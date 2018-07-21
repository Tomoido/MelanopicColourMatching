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
wls_snb = mat(:, 1);
snb_rgb = mat(:, 2:4);
r_match = round(645.16);
g_match = round(526.32);
b_match = round(444.44);
primaries = [r_match,g_match,b_match];

% load cone sensitivities
wls = (390:1:830)'; 
S = WlsToS(wls);
receptorObjFieldSize10Deg = SSTReceptorHuman('S', S, 'fieldSizeDeg', 10);
T_receptors = receptorObjFieldSize10Deg.T.T_energyNormalized; % starting with normlized
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);

% set up parameters for minimization algorithms
n_test = length(wls_snb);
startp = [-1,1,1]; 
minp = [-5,-5,-5]; 
maxp = [5,5,5];
options = optimoptions('fmincon','Display', 'none',...
    'OptimalityTolerance',      0,...
    'StepTolerance',            0,...
    'MaxIterations',            3e3,...
    'MaxFunctionEvaluations',   3e3);


% algorithm minimizing squared difference in normalized stimulation
for i=1:n_test
    [ndiff_rgb(i,:), ndiff_fval(i)] = fmincon(@(T_rgb)opt_primaries_diff...
        (T_rgb,wls(i),primaries,[L;M;S]),startp,...
        [],[],[],[],minp,maxp,[], options);
end


% algorithm minimizing squared difference in absolute stimulation
% loading non-normalized cone sensitivities
T_receptors = receptorObjFieldSize10Deg.T.T_energy;
L = T_receptors(1,:);
M = T_receptors(2,:);
S = T_receptors(3,:);
mel = T_receptors(4,:);

for i=1:n_test
    [diff_rgb(i,:), diff_fval(i)] = fmincon(@(T_rgb)opt_primaries_diff...
        (T_rgb,wls(i),primaries,[L;M;S]),startp,...
        [],[],[],[],minp,maxp,[], options);
end


% algorithm minimizing squared cone contrast
% adding more starting values: where each primary is negative once
startp = [-1,1,1; 1,-1,1; 1,1,-1];
for i=1:n_test
    for c = 1:3 % we go through this extra loop to load the extra starting values
        % this ensures that the algorithm won't get stuck in local minima
        % and will converge in a reasonable number of
        % iterations,f evaluations
        [return_rgb(i,:,c),fval(i,c)] = fmincon(@(T_rgb)opt_primaries_cont...
            (T_rgb,wls(i),primaries,[L;M;S]),startp(c,:),...
            [],[],[],[],minp,maxp,[], options);
    end
    [cont_fval(i),position(i)] = min(fval(i,:));
    cont_rgb(i,:) = return_rgb(i,:,position(i));
end

% saving output CMFs
save('algorithm_CMFs','cont_rgb','diff_rgb','ndiff_rgb')


%% Plot output CMFs against SNB data

figure

subplot(3,1,1)
title('using sum of squared cone contrast')
r = plot(wls_snb,snb_rgb(:,1),'r','LineWidth',5); hold on
rm = plot(wls_snb,cont_rgb(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(wls_snb,snb_rgb(:,2),'g','LineWidth',5); hold on
gm = plot(wls_snb,cont_rgb(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(wls_snb,snb_rgb(:,3),'b','LineWidth',5); hold on
bm = plot(wls_snb,cont_rgb(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;
legend({'r SNB','r model','g SNB','g model','b SNB','b model'})

subplot(3,1,2)
title('using sum of squared differences (non-normalized)')
r = plot(wls_snb,snb_rgb(:,1),'r','LineWidth',5); hold on
rm = plot(wls_snb,diff_rgb(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(wls_snb,snb_rgb(:,2),'g','LineWidth',5); hold on
gm = plot(wls_snb,diff_rgb(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(wls_snb,snb_rgb(:,3),'b','LineWidth',5); hold on
bm = plot(wls_snb,diff_rgb(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;

title('using sum of squared differences (normalized)')
subplot(3,1,3)
r = plot(wls_snb,snb_rgb(:,1),'r','LineWidth',5); hold on
rm = plot(wls_snb,ndiff_rgb(:,1),'r:','LineWidth',2);
r.Color(4) = 0.25;
rm.Color(4) = 0.75;
g = plot(wls_snb,snb_rgb(:,2),'g','LineWidth',5); hold on
gm = plot(wls_snb,ndiff_rgb(:,2),'g:','LineWidth',2);
g.Color(4) = 0.25;
gm.Color(4) = 0.75;
b = plot(wls_snb,snb_rgb(:,3),'b','LineWidth',5); hold on
bm = plot(wls_snb,ndiff_rgb(:,3),'b:','LineWidth',2);
b.Color(4) = 0.25;
bm.Color(4) = 0.75;


%% Function for the difference minimization algorithms

function sqdiff = opt_primaries_diff(T_rgb,wls_test,primaries,cone_sens)

if sum(T_rgb<0)~=1
    sqdiff = inf;   % catch if more than one primary is negative
else
    
    % sort negative and positive primaries
    neg_prim = find(T_rgb<0);
    pos_prim(:) = find(T_rgb>0);
    
    wls = (390:1:830)';  % wavelengths for spectra
    n_wls = length(wls);
    spdt = zeros(n_wls,1);  % test spectrum var
    spdm = zeros(n_wls,1);  % match spectrum var

    FWHM = 10; % half width of gaussian based on interference filter properties
    sigma = FWHM/2.4; 
    gauss = @(height,position)height.*exp(-((wls-position).^2)/(2*sigma^2));
      
    % construct test spectrum
    spdt = spdt+gauss(1,wls_test);
    spdt = spdt+gauss(abs(T_rgb(neg_prim)),primaries(neg_prim));

    % construct match spectrum
    spdm = spdm+gauss(T_rgb(pos_prim(1)),primaries(pos_prim(1)));
    spdm = spdm+gauss(T_rgb(pos_prim(2)),primaries(pos_prim(2)));
    
    % calculate cone stimulation
    stimm = cone_sens(1:3,:)*spdm;
    stimt = cone_sens(1:3,:)*spdt;

    % calculate sum of squared differences in stimulation
    sqdiff = sum((stimt-stimm).^2);
    
end
end

%% Function for the contrast minimization algorithms

function sqcont = opt_primaries_cont(T_rgb,wls_test,primaries,cone_sens)

if sum(T_rgb<0)~=1
    sqcont = inf;   % catch if more than one primary is negative    
else
    
    % sort negative and positive primaries
    neg_prim = find(T_rgb<0);
    pos_prim(:) = find(T_rgb>0);
    
    wls = (390:1:830)';  % wavelengths for spectra
    n_wls = length(wls);
    spdt = zeros(n_wls,1);  % test spectrum var
    spdm = zeros(n_wls,1);  % match spectrum var

    FWHM = 10; % half width of gaussian based on interference filter properties
    sigma = FWHM/2.4; 
    gauss = @(height,position)height.*exp(-((wls-position).^2)/(2*sigma^2));
      
    % construct test spectrum
    spdt = spdt+gauss(1,wls_test);
    spdt = spdt+gauss(abs(T_rgb(neg_prim)),primaries(neg_prim));

    % construct match spectrum
    spdm = spdm+gauss(T_rgb(pos_prim(1)),primaries(pos_prim(1)));
    spdm = spdm+gauss(T_rgb(pos_prim(2)),primaries(pos_prim(2)));
    
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
