function [seg] = getPermeability(MP_Constant,seg)
%getPermeability calculates the fraction of the input signal that is lost
%as the signal passes through the vessel, based on the vessel length and
%the permeability of the vessel per unit length
%   INPUTS: seg structure produced by sortSegments
%   OUTPUTS: seg structure with appended field, permeability

%to convert form cm/s to um/s
conversionFactor = 1e4;

for i= 1:numel(seg)
    seg(i).MP = MP_Constant*conversionFactor*seg(i).Length*2*pi*seg(i).Radius*seg(i).Delay/seg(i).Volume;
    
    if seg(i).MP > 1
        seg(i).MP = 1;
    end
end

end