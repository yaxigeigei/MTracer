classdef regressionMAELayer < nnet.layer.RegressionLayer
               
    methods
        function layer = regressionMAELayer(name)
            % Create an exampleRegressionMAELayer

            % Set layer name
            if nargin == 1
                layer.Name = name;
            end

            % Set layer description
            layer.Description = 'Regression layer with MAE loss';
        end
        
        function loss = forwardLoss(layer, Y, T)
            % Returns the MAE loss between the predictions Y and the
            % training targets T

            % Calculate MAE
            K = size(Y,3);
            meanAbsoluteError = sum(abs(Y-T),3)/K;
    
            % Take mean over mini-batch
            N = size(Y,4);
            loss = sum(meanAbsoluteError)/N;
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            % Returns the derivatives of the MAE loss with respect to the predictions Y

            N = size(Y,4);
            R = size(Y,3);
            dLdY = sign(Y-T)/N/R;
        end
    end
end