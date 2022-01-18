classdef ContourInfo < handle
properties
    total_cnt = 0;
    contour_store = {};
    llr_contour_store = {};
    llr_interp_store = {};
    weight_component_store = {};
end

methods
    function obj = ContourInfo()
    end
    
    function add_contour(obj, rot_contour, weight_cmp, llr_interp)
        obj.total_cnt = obj.total_cnt + 1;
        obj.contour_store{obj.total_cnt} = rot_contour;
        obj.llr_interp_store{obj.total_cnt} = llr_interp;
        obj.weight_component_store{obj.total_cnt} = weight_cmp;
        if size(rot_contour, 2) == 4
            obj.llr_contour_store{obj.total_cnt} = geo.quat2llr(rot_contour);
        end
    end
    
    display_info(obj)  % Defined in a seperated file.
end
end