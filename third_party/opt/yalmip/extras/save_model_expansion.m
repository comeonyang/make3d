function done = save_model_expansion(method,F_graph,properties)
global ALREADY_MODELLED REMOVE_THESE_IN_THE_END

done = 0;

for variables = properties.models
    if length(ALREADY_MODELLED)>=variables & ~isempty(ALREADY_MODELLED{variables})
        % Ouch, we have already modelled this-one
        if strcmpi(ALREADY_MODELLED{variables}.method,method)
            % Ok, already modelled using same approach
            done = 1;
            return
        elseif strcmpi(ALREADY_MODELLED{variables}.method,'graph') & strcmpi(method,'milp')
            % Replace old graph model with MILP model
            REMOVE_THESE_IN_THE_END = [REMOVE_THESE_IN_THE_END  ALREADY_MODELLED{variables}.index];
            ALREADY_MODELLED{variables}.method = 'milp';
            ALREADY_MODELLED{variables}.index  = getlmiid(F_graph);
        elseif strcmpi(ALREADY_MODELLED{variables}.method,'milp') & strcmpi(method,'graph')
            % Keep old stuff, we are done
            done = 1;
            return
        end
    else
        ALREADY_MODELLED{variables}.method = method;
        ALREADY_MODELLED{variables}.index  = getlmiid(F_graph);
    end
end
