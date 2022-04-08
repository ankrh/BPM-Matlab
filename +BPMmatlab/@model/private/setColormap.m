function setColormap(hAxis,colormapType)
switch colormapType
  case 'GPBGYR'
    colormap(hAxis,GPBGYRcolormap);
  case 'HSV'
    colormap(hAxis,hsv/1.5);
  case 'Parula'
    colormap(hAxis,parula);
  case 'Gray'
    colormap(hAxis,gray);
  case 'Cividis'
    colormap(hAxis,cividisColormap);
end
end