function x = perturbPt(x, bounds, range)
  numDims = size(bounds, 1);
  perturbation =  -range + 2 * range * rand(1, numDims);
  x = x + perturbation;
  x = projectToRectangle(x', bounds)';
end
