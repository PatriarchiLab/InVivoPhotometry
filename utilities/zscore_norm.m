function [out]= zscore_norm(in, range1, range2)
    mu=mean(in(range1:range2));
    sigma=std(in(range1:range2));
    out=(in-mu)/sigma;
end

