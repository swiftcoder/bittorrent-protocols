% =========================================================================
% Purpose: This M-File Script creates theorical exponential distribution
%          for the P2P Protocols.
%
% Support: exponential(lambda, n)
%           exponential_distribution(lambda, n, filename)          
%
% Date   : 03/05/10
% Author : Tommy Nguyen
% =========================================================================

% Return a matrix containing n values and saves values in a text file.
function exponential = exponential_distribution(lambda, n, filename)
    fid = fopen(filename, 'w');
    for j=1:n,
        exponential_value = -log(rand(1))./(1/lambda);
        fprintf(fid, '%f \n', exponential_value);
        exponential{j} = exponential_value;
    end
end
