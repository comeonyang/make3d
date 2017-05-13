% RANSAC - Robustly fits a model to data with the RANSAC algorithm
%
% Usage:
%
% [M, inliers] = ransac(x, fittingfn, distfn, degenfn s, t, feedback)
% [F, inliers, NewDistrib, fail] = ...
% ransac(defaultPara, x, Depth, ...
%        fittingfn, distfn, degenfn, s, t, distrib, T, feedback, disp, FlagDist);
%
% Arguments:
%     defaultPara - useful default parameters (like, camera intrinsic matrix)
%
%     x         - Data sets to which we are seeking to fit a model M
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
%
%     Depth     - depth imformation to support more accurate ransac
%
%     fittingfn - Handle to a function that fits a model to s
%                 data from x.  It is assumed that the function is of the
%                 form: 
%                    M = fittingfn(x)
%                 Note it is possible that the fitting function can return
%                 multiple models (for example up to 3 fundamental matrices
%                 can be fitted to 7 matched points).  In this case it is
%                 assumed that the fitting function returns a cell array of
%                 models.
%                 If this function cannot fit a model it should return M as
%                 an empty matrix.
%
%     distfn    - Handle to a function that evaluates the
%                 distances from the model to data x.
%                 It is assumed that the function is of the form:
%                    [inliers, M] = distfn(M, x, t)
%                 This function must evaluate the distances between points
%                 and the model returning the indices of elements in x that
%                 are inliers, that is, the points that are within distance
%                 't' of the model.  Additionally, if M is a cell array of
%                 possible models 'distfn' will return the model that has the
%                 most inliers.  If there is only one model this function
%                 must still copy the model to the output.  After this call M
%                 will be a non-cell object representing only one model. 
%
%     degenfn   - Handle to a function that determines whether a
%                 set of datapoints will produce a degenerate model.
%                 This is used to discard random samples that do not
%                 result in useful models.
%                 It is assumed that degenfn is a boolean function of
%                 the form: 
%                    r = degenfn(x)
%                 It may be that you cannot devise a test for degeneracy in
%                 which case you should write a dummy function that always
%                 returns a value of 1 (true) and rely on 'fittingfn' to return
%                 an empty model should the data set be degenerate.
%
%     s         - The minimum number of samples from x required by
%                 fittingfn to fit a model.
%
%     t         - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%
%     distrib  - initial distribution (default uniform dist)
%
%     T - affine transform 3 by 3 matrix applied on x
%
%     feedback  - An optional flag 0/1. If set to one the trial count and the
%                 estimated total number of trials required is printed out at
%                 each step.  Defaults to 0.
%
%     disp  - if true, display the matches found when done.
%
%     FlagDist - if true, calculate the reprojection error
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
%     NewDist - New Distribution after Ransac (Outliers have zero distribution)
%     fail    - true if Ransac fail to find any solution is not degenerated
%
% For an example of the use of this function see RANSACFITHOMOGRAPHY or
% RANSACFITPLANE 

% References:
%    M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
%    for model fitting with applications to image analysis and automated
%    cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981
%
%    Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
%    Computer Vision". pp 101-113. Cambridge University Press, 2001

% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au    
% http://www.csse.uwa.edu.au/~pk
%
% May      2003 - Original version
% February 2004 - Tidied up.
% August   2005 - Specification of distfn changed to allow model fitter to
%                 return multiple models from which the best must be selected.
% 
% additional parameter distrib is a vector of non-negative numbers that
% specifies a (not necessarily normalized) probability distribution over
% the different possible matches - passed to the sample function in the lightspeed
% package in order to sample possible matches. 
% (added by Jeff Michels)
% Additional NewDist estimated after Ransac is formed by using Depth
% information
% (added by Min Sun)
%function [M, inliers, fail] = ransac(x, fittingfn, distfn, degenfn, s, t, feedback, distrib)
function [M, inliers, NewDistrib, fail] = ...
          ransac(defaultPara, x, Depth, ...
          fittingfn, distfn, degenfn, s, t, distrib, T, feedback, disp, FlagDist)
    
    if nargin < 10
	feedback = 0;
        disp = 0;
    end
    fail=0;
    [rows, npts] = size(x);                 
    
    %p = 0.999;    % Desired probability of choosing at least one sample
    p = 0.9999;    % Desired probability of choosing at least one sample
                 % free from outliers

    %maxTrials = 20000;    % Maximum number of trials before we give up.
    %maxDataTrials = 100; % Max number of attempts to select a non-degenerate
    maxTrials = 20000;    % Maximum number of trials before we give up.
    maxDataTrials = 100; % Max number of attempts to select a non-degenerate
                         % data set.
    
    bestM = NaN;      % Sentinel value allowing detection of solution failure.
    trialcount = 0;
    bestscore =  0;    
    N = 1;            % Dummy initialisation for number of trials.
    
    while N > trialcount
        
        % Select at random s datapoints to form a trial model, M.
        % In selecting these points we have to check that they are not in
        % a degenerate configuration.
        degenerate = 1;
        count = 1;
        while degenerate
            % Generate s random indicies in the range 1..npts
            %ind = ceil(rand(1,s)*npts);
            % JM replaced above line with sample() to sample non-uniformly
            ind = sample(distrib, s);
            
            % Test that these points are not a degenerate configuration.
            degenerate = feval(degenfn, x(:,ind));
	    
	    if ~degenerate 
		% Fit model to this random selection of data points.
		% Note that M may represent a set of models that fit the data in
		% this case M will be a cell array of models
		M = feval(fittingfn, x(:,ind));
		
		% Depending on your problem it might be that the only way you
                % can determine whether a data set is degenerate or not is to
                % try to fit a model and see if it succeeds.  If it fails we
                % reset degenerate to true.
		if isempty(M)
		    degenerate = 1;
		end
	    end
	    
	    % Safeguard against being stuck in this loop forever
	    count = count + 1;
            if count > maxDataTrials
                warning('Unable to select a nondegenerate data set');
                break
            end
        end
        
        % Once we are out here we should have some kind of model...        
        % Evaluate distances between points and model returning the indices
        % of elements in x that are inliers.  Additionally, if M is a cell
        % array of possible models 'distfn' will return the model that has
        % the most inliers.  After this call M will be a non-cell object
        % representing only one model.
        [inliers, M, tempReProjError] = ...
                  feval(distfn, M, x, t, defaultPara, Depth, T, FlagDist);
	warning('called distfn in ransac');
        if ~isempty(tempReProjError)
            ReProjErrorM(trialcount+1,:) = tempReProjError;
        else
            ReProjErrorM = [];
        end
        
        % Find the number of inliers to this model.
        ninliers = length(inliers);
        
        % weightedNInliers weights each inlier by the prior probability
        % that it should be an inlier based on the parameter distrib
        % added by JM
        weightedNInliers = sum(distrib(inliers));
        
        %if ninliers > bestscore    % Largest set of inliers so far...
        %    bestscore = ninliers;  % Record data for this model
        if weightedNInliers > bestscore
            bestscore = weightedNInliers;
            bestinliers = inliers;
            bestM = M;
            
            % Update estimate of N, the number of trials to ensure we pick, 
            % with probability p, a data set with no outliers.
            fracinliers =  weightedNInliers;%ninliers/npts;
            pNoOutliers = 1 -  fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
        end
        
        trialcount = trialcount+1;
	if feedback&&mod(trialcount,100)==0&& trialcount>1
	    fprintf('trial %d out of %d         \r',trialcount, ceil(N));
	end

        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            warning( ...
            sprintf('ransac reached the maximum number of %d trials',...
                    maxTrials));
            fail=1;
            break
        end     
    end
    fprintf('\n');
    
    if ~isnan(bestM)   % We got a solution 
        M = bestM;
        inliers = bestinliers;
    else           
        warning('ransac was unable to find a useful solution');
        fail=1;
    end
   
    % Processing NewDist
    if FlagDist
        ReProjError = median(ReProjErrorM,1);
        if disp
            figure; hist(ReProjError(ReProjError ~= Inf)); % plot hist of the ReProjection error
        end
        VarGammDist =  var(ReProjError(ReProjError ~= Inf));
        NewDistrib(size(ReProjError)) = 0;
        NewDistrib(ReProjError(ReProjError ~= Inf)) = ...
                   gampdf(ReProjError(ReProjError ~= Inf), 1, VarGammDist);
    else
        NewDistrib = [];
    end
