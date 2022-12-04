function stop = save_pso(optimValues, state, prefix) 
    save([prefix, num2str(optimValues.iteration),'.mat'], 'optimValues','-v7.3');
    stop = false;
end
