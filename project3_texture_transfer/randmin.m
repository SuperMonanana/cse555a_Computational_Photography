function [val, ind] = randmin(x);
  
  perm = randperm(length(x));
  [val ind] = min(x(perm));
  ind = perm(ind);