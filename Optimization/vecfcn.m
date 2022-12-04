function y = vecfcn(pop, fun, pool, t_threshold)
    clear future;
    
    popSize = size(pop,1);
    y = nan(popSize,1);
    future_index = zeros(2*pool.NumWorkers,1);
    
    nbytes = fprintf('Finished 0 of %d particles', popSize);
    ntimeout = 0;
    
    for i = 1:min(2*pool.NumWorkers,popSize) % your outer for-loop()
        % Invoke your function on a worker
        future_index(i) = i; 
        future(i) = parfeval(pool, fun, 1, pop(i,:));
    end
    
    done=repmat(false,1,popSize);
    while ~all(done)
        fprintf(repmat('\b',1, nbytes));
        nbytes = fprintf('Finished %d of %d particles (%d timeouts)', sum(done), popSize, ntimeout);
        
        for id = 1:size(future,2)
            f = future(id);
            initialize=false;
            
            if ~isempty(f.Error)
                % disp(f);
                % disp(pop(future_index(id),:));
                y(future_index(id),:) = NaN;
                initialize = true; 
            elseif strcmp(f.State, 'finished') & ~f.Read & isempty(f.Error)
                y(future_index(id),:) = fetchOutputs(f);
                initialize=true;
            elseif strcmp(f.State, 'running') 
                if datetime(datetime, 'TimeZone', 'local') - f.StartDateTime >= t_threshold
                    ntimeout = ntimeout + 1;
                    cancel(f);
                    y(future_index(id),:) = NaN;
                    initialize = true;
                end
            end
            
            if initialize
                done(future_index(id)) = true; 
                i = i + 1;
                if i <= popSize
                    future_index(id) = i; 
                    future(id) = parfeval(pool, fun, 1, pop(i,:));
                end
            end
        end
    end
    
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('Finished %d of %d particles (%d timeouts)\n', sum(done), popSize, ntimeout);

    clear future;
end
