% =========================================================================
% Purpose: This M-File Script creates a theorical log-normal distribution 
%          P2P Protocols. 
%
% Support: lognorma_distribution(mu, sigma, n, filename)
%                   
% Date   : 03/06/10
% Author : Tommy Nguyen
% =========================================================================

% Function returns a matrix containing n values that follows a log-normal
% distribution. The data is also saved into a file. 

function lognormal = lognormal_distribution(mu, sigma, n, filename)
    fid = fopen(filename, 'w');
    for i=1:n,
        lognormal_value = lognrnd(mu,sigma)
        fprintf(fid, '%f \n', lognormal_value);
        lognormal{i} = lognormal_value;
    end
end