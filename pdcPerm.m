function [bcpR, bcpP, permR, permP ,rid] = pdcPerm(x,y,z,perms,para)
% This function interfaces with pdc.m to produce the bias corrected partial
% distance correlation coefficient as well as a permutation-based p-value.
% A single variable is permuted (y) a number of times (perms) and the pdc
% is computed. The proportion of permutations where permuted pdc exceeds
% actual pdc is taken as the permuted p-value.

% INPUTS:
% x and y: n x 1 vectors that you are interested in understanding the
%   relationship between
% z: n x p matrix of p variables to condition the relationship between x
%   and y.
% perms: number of permutations to perform (recommended: 10000, or at least
%   1000). This can be set to 0 to perform no permutations.
% para: if true, will use parfor for parellalization.
%
% OUTPUTS:
% bcpR: bias corrected partial distance correlation "coefficient"
% bcpP: bias corrected partial distance correlation p-value (application of
%   standard bias corrected distance correlation method to partial distance correlation; see pdc.m)
% permR: permuted bcpR values
% permP: permutation-based p-value for bcpR
% rid: samples that were removed because they were identified to be NaNs
%
% Example calls : 
% [bcpR, bcpP, permR, permP, rid] = pdcPerm(x,y,z,10,true);
% [bcpR, bcpP, permR, permP, rid] = pdcPerm(x,y,z,0,true);
% [bcpR, bcpP, permR, permP, rid] = pdcPerm(x,y,z,[],[]);
%
% alex.teghipco@sc.edu

if isempty(x) || isempty(y) || isempty(z)
    error('x, y or z variable inputs are empty...please provide data')
end

if isempty(perms)
    perms = 10000;
end

if isempty(para)
    para = true;
    warning('Running parallelized for loop by default, if pdcPerm crashes please set input variable para to false...')
end

if min(size(x)) > 1 || min(size(y)) > 1
    error('x or y variables have too many columns...samples are rows and both x and y should contain one variable/column...')
end

if size(x,2) > size(x,1)
    x=x';
    warning('Flipping x matrix so that rows are columns and columns are rows, double check that your inputs are correct and n x 1 with n samples')
end

if size(y,2) > size(y,1)
    y=y';
    warning('Flipping y matrix so that rows are columns and columns are rows, double check that your inputs are correct and n x 1 with n samples')
end

if size(z,2) > size(z,1)
    warning('It looks like the z matrix is incorrect but I canot be certain so please check yourself. z should be an n x p matrix with n (rows) samples and p (cols) control variables')
end

id1 = find(isnan(x));
id2 = find(isnan(y));
if size(z,2) > 1
    [id3, ~] = find(isnan(z));
else
    id3 = find(isnan(z));
end
rid = unique([id1; id2; id3]);
x(rid) = [];
y(rid) = [];
if size(z,2) > 1
    z(rid,:) = [];
else
    z(rid) = [];
end
if ~isempty(rid)
    warning(['removed ' num2str(length(rid)) ' samples...check output variable rid'])
end

[bcpR, bcpP] = pdc(x,y,z,'distance',true);
if perms > 0
    rng('shuffle');
    if para
        ppm = waitbar(0,'Please wait ...');
        p = 1;
        D = parallel.pool.DataQueue;
        afterEach(D, @nUpdateWaitbar);
        parfor pi = 1:perms
            yperm = y(randperm(size(y,1)));
            [permR(pi,1), ~] = pdc(x,yperm,z,'distance',true);
            send(D,p)
        end
        close(ppm)
    else
        for pi = 1:perms
            disp(['Working on permutation ' num2str(pi) ' of ' num2str(perms)])
            yperm = y(randperm(size(y,1)));
            [permR(pi,1), ~] = pdc(x,yperm,z,'distance',true);
        end
    end
    permP = length(find(permR >= bcpR))/perms;
    

else
    permR = [];
    permP = [];   
end

    function nUpdateWaitbar(~)
    waitbar(p/perms, ppm,['Working on permutation ' num2str(p) ' of ' num2str(perms)]);
    p = p + 1;
    end
end