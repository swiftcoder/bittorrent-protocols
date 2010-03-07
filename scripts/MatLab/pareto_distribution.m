% =========================================================================
% Purpose: This M-File Script creates theorical pareto distribution
%          for the P2P Protocols.
%
% Support: pareto(alpha, n, filename)
%                   
% Date   : 03/05/10
% Author : Tommy Nguyen
% =========================================================================

% Function returns a matrix containing n values that follows a pareto
% distribution. The data is also saved into a file. 

function pareto = pareto_distribution(alpha, n, filename)
    fid = fopen(filename, 'w');
    for i=1:n,
        pareto_value = (1 - rand(1)).^(-1/alpha)-1;
        fprintf(fid, '%f \n', pareto_value);
        pareto{i} = pareto_value;
    end
end
