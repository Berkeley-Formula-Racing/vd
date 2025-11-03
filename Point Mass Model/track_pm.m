
classdef track_pm
    properties
        s double  % arclength [m], size N
        k double  % curvature [1/m], size N
    end
    methods
        function obj = track_pm(s, k)
            obj.s = s(:);
            obj.k = k(:);
        end
    end
    methods (Static)
        function trk = loadFromMichiganMat(matfile)
            data = load(matfile);
            s = data.arclength(:);
            k = data.curvature(:);
            trk = track_pm(s, k);
        end
    end
end
