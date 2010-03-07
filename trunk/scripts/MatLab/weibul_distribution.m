% =========================================================================
% Purpose: This M-File Script creates a theoretical Weibul Distribution. 
%
% Support: weibul_distribution(scale, shape, n, file)
%                   
% Date   : 03/05/10
% =========================================================================

% Function returns a matrix contain n values of the Weibul Distribution

function weibul = weibul_distribution(scale, shape, n, filename)
    fid = fopen(filename, 'w');
    for i=1:n,
        weibul_value = scale.*(-log(1-rand(1,1))).^(1/shape);
        fprintf(fid, '%f \n', weibul_value);
        weibul{i} = weibul_value;
    end
end