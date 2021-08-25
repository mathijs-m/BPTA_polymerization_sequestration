function [stop] = interuptFun(x,optimValues,state,varargin,best_res,max_iter)
%Interuption function to check if lsqnonlin achieves a decently low
%residual after max_iter iterations. If the residual is not in the same
%ball park as the beste achieved one, then stop the optimization, since it
%is going nowhere, so don't waste time.
if isequal(state,'init') && optimValues.resnorm >=400000
    stop = 1;
elseif isequal(state,'iter') && optimValues.iteration >= 5 && optimValues.resnorm >=36000
    stop = 1;
elseif isequal(state,'iter') && optimValues.iteration >= 15 && optimValues.resnorm >=100000
    stop = 1;
elseif isequal(state,'iter') && optimValues.iteration >= max_iter && optimValues.resnorm >=1.3*best_res
    stop = 1;
else
    stop = 0;
end

end

