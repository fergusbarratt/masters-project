function [] = iphnum_vs_purity(EVals, DVals, varargin)
  %Must close gcf beforehand
  if nargin>2
    both = varargin{1};
  else
    both = 0;
  end
  % delete(gcf)
  if both
    subplot(1, 2, 1)
    prettyplot(EVals, DVals, iphnumvals(EVals, DVals, 'iphnum'), 'Iphnum');
    subplot(1, 2, 2)
    prettyplot(EVals, DVals, iphnumvals(EVals, DVals, 'purity'), 'Purity');
  else
    prettyplot(EVals, DVals, iphnumvals(EVals, DVals, 'iphnum', 60), 'Iphnum');
  end
