% CONE STIMULATION DIFFERENCE FUNCTION FOR CMF GENERATION ALGORITHM
%------------------------------------------------------------------

% This function returns the sum of squared differences in stimulation of
% the three cone types. It is necessary for the CMF generation function.

function sqdiff = GenerateCMF(CMF,wls_test,primaries,cone_sens,range)

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
