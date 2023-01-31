function [vec_out] = denanvec(vec_in)

vec_out = vec_in(~isnan(vec_in));

end