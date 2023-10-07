%------Functions -------%
function [err] = loss_function(w, w_experimental, type) % or cost_function
   c = Constants; % Gets constants for the model
   % First normalise the objective function
   d = w - w_experimental;
   %d = d / max(abs(d));
    % Log-Cosh Loss
    % https://heartbeat.fritz.ai/5-regression-loss-functions-all-machine-learners-should-know-4fb140e9d4b0
    % https://memoex.github.io/note/tech/ml/loss/
    if type == 1
        err = 0.0;
        for i=1:1:length(d)
            err = err + log(cosh(d(i)));
        end
        if c.loss_function_division == true
            err = err / length(d);
        end
        return;
    end
    % MSE (Mean square error)
    if type == 2
        err = norm(d)^2;
        if c.loss_function_division == true
            err = err / length(d);
        end
        return;
    end
    
    % RMSE (Root mean square error)
    % https://www.vernier.com/til/1014
    % https://medium.com/human-in-a-machine-world/mae-and-rmse-which-metric-is-better-e60ac3bde13d
    if type == 3
        err = norm(d);
        if c.loss_function_division == true
            err = err / sqrt(length(d));
			%err = err / length(d);
        end
        return;
    end
end
%--------------------%