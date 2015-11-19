function [plo] = prettyplot(Xrange, Yrange, h, varargin)
  if isrow(h)
    plot(Xrange, h)
  elseif iscolumn(h)
    plot(Yrange, h)
  else 
    plo = surf(Xrange, Yrange, h);
    plo.LineStyle = 'none';
    view(3);
    colormap('jet');
    colorbar
    % ylim=get(gca,'YLim');
    % xlim=get(gca,'XLim');
    % text((xlim(1)-1),(ylim(1)-1),[num2str(E) ',' num2str(det)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
    % xlabel('Re(Q)')
    % ylabel('Im(Q)')
  end
    if nargin>3
      try
        title(varargin{1})
      catch
        'title not a string'
      end
    end
