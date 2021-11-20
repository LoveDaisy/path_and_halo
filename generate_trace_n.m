function trace_n = generate_trace_n(crystal, trace)
trace_n = crystal.n * ones(size(trace.fid));
for j = 1:length(trace_n)
    trace_n(j) = trace_n(j) * (-1)^(j+1);
end
if length(trace_n) > 1
    trace_n(end) = sign(trace_n(end - 1));
end
trace_n = [1; trace_n];
end