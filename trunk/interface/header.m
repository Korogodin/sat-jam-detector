if Kind == 1 
    hF = figure(hObject);
    hObject = axes('Parent', hF);
else
    BDF = get(hObject, 'ButtonDownFcn');
end
hA = hObject;