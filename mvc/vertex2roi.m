function roi = vertex2roi( atlas, vertex )
scouts = atlas.Scouts;
roi = NaN;
for vi = 1 : numel(scouts)
    scoutsearch = scouts(vi).Vertices;
    if any(scoutsearch == vertex)
        roi = vi;
    end
end
end

