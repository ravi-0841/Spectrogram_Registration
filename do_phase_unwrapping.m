%% This function unwraps the phase to give smoother phase plot.
function unwrapped_phase = do_phase_unwrapping(phase)
    unwrapped_phase = phase;
    for col = 1:size(phase,2)
        for row = 2:size(phase,1)
            if (unwrapped_phase(row,col) - unwrapped_phase(row-1,col)) > pi
                unwrapped_phase(row:end,col) = unwrapped_phase(row:end,col) - 2*pi;
            elseif (unwrapped_phase(row,col) - unwrapped_phase(row-1,col)) < -1*pi
                unwrapped_phase(row:end,col) = unwrapped_phase(row:end,col) + 2*pi;
            end
        end
    end
end