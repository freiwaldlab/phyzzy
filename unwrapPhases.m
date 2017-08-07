function [ phases, piMultiples ] = unwrapPhases( phases )
%unwrapPhases lets you plot intuitive phase plots without polar axes
%   by undoing the mod pi operation in the most intuitive way

piMultipliers = zeros(size(phases));
for spec_i = 1:size(phases,1) %spec short for spectrum
  for freq_i = 2:size(phases,2)
    if phases(spec_i,freq_i-1) < -pi/2 && phases(spec_i,freq_i) > pi/2
      piMultipliers(spec_i,freq_i) = piMultipliers(spec_i,freq_i-1) - 1;
      continue
    end
    if phases(spec_i,freq_i-1) > pi/2 && phases(spec_i,freq_i) < -pi/2
      piMultipliers(spec_i,freq_i) = piMultipliers(spec_i,freq_i-1) + 1;
      continue
    end
    piMultipliers(spec_i,freq_i) = piMultipliers(spec_i,freq_i-1);
  end
end
phases = phases + 2*pi*piMultipliers;
piMultiples = ceil(min(min(phases))/pi):floor(max(max(phases))/pi);
end

