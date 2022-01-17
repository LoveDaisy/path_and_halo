classdef Progress < handle
properties (Access = private)
    curr_prg = 0;
    step = 0.01;
    update_prg = 0;
    update_pct_step = 0.05;
end

methods
    function obj = Progress(step, update_pct_step)
    if nargin == 0
    elseif nargin == 1
        obj.step = step;
    elseif nargin == 2
        obj.step = step;
        obj.update_pct_step = update_pct_step;
    else
        error('Too many input arguments!');
    end
    end

    function tik(obj, n_steps)
    obj.curr_prg = obj.curr_prg + n_steps * obj.step;
    obj.update_prg = obj.update_prg + n_steps * obj.step;
    end

    function show(obj, msg_str, varargin)
    fprintf('[%05.2f%%]', obj.curr_prg * 100);
    if ~isempty(varargin)
        fprintf(msg_str, varargin{:});
    end
    fprintf('\n');
    end

    function res = check_update(obj)
    res = obj.update_prg > obj.update_pct_step;
    obj.update_prg = mod(obj.update_prg, obj.update_pct_step);
    end

    function need_update = tik_and_show(obj, n_steps, msg_str, varargin)
    obj.tik(n_steps);
    need_update = obj.check_update();
    if need_update
        obj.show(msg_str, varargin{:});
    end
    end
end
end
