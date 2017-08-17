function [final_vps, final_NFAs] = remove_duplicates3(vps, NFAs, params)
% identifies and removes duplicate detections, keeping only most
% significant ones.

THRESHOLD = params.DUPLICATES_THRESHOLD;


clus = aggclus(vps', THRESHOLD);


final_vps = [];
final_NFAs = [];
for i=1:length(clus)
    c = clus{i};
    if length(c)==1
        final_vps = [final_vps, vps(:,c)];
        final_NFAs = [final_NFAs; NFAs(c)];
    else
        [V ,I] = min(NFAs(c));
        final_vps = [final_vps, vps(:,c(I))];
        final_NFAs = [final_NFAs; NFAs(c(I))];
    end
    
end

