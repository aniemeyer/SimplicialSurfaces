#############################################################################
##
##  SimplicialSurface package
##
##  Copyright 2012-2019
##    Markus Baumeister, RWTH Aachen University
##    Alice Niemeyer, RWTH Aachen University 
##
## Licensed under the GPL 3 or later.
##
#############################################################################


#######################################
##
##      Moving along edges
##
InstallMethod( IsVerticesAdjacentNC, "for a VEF-complex and two vertices",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, v1, v2)
        if v1 < v2 then
            return [v1,v2] in VerticesOfEdges(complex);
        else
            return [v2,v1] in VerticesOfEdges(complex);
        fi;
    end
);
InstallMethod( IsVerticesAdjacent, "for a VEF-complex and two vertices",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, v1, v2)
        local name;

        name := "IsVerticesAdjacent";
        __SIMPLICIAL_CheckVertex(complex, v1, name);
        __SIMPLICIAL_CheckVertex(complex, v2, name);

        return IsVerticesAdjacentNC(complex, v1, v2);
    end
);

InstallMethod( EdgesBetweenVerticesNC, "for a VEF-complex and two vertices",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, v1, v2)
        local res, e1, e2, e;

        res := [];
        e1 := EdgesOfVertices(complex)[v1];
        e2 := EdgesOfVertices(complex)[v2];
        for e in e1 do
            if e in e2 then
                Add(res,e);
            fi;
        od;

        return res;
    end
);
InstallMethod( EdgesBetweenVertices, "for a VEF-complex and two vertices",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, v1, v2)
        local name;

        name := "EdgesBetweenVertices";
        __SIMPLICIAL_CheckVertex(complex, v1, name);
        __SIMPLICIAL_CheckVertex(complex, v2, name);

        return EdgesBetweenVerticesNC(complex, v1, v2);
    end
);

InstallMethod( EdgeBetweenVerticesNC, "for a VEF-complex and two vertices",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, v1, v2)
        local edges;

        edges := EdgesBetweenVerticesNC(complex, v1, v2);
        if Length(edges) = 1 then
            return edges[1];
        else
            return fail;
        fi;
    end
);
InstallMethod( EdgeBetweenVertices, "for a VEF-complex and two vertices",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, v1, v2)
        local name;

        name := "EdgeBetweenVertices";
        __SIMPLICIAL_CheckVertex(complex, v1, name);
        __SIMPLICIAL_CheckVertex(complex, v2, name);

        return EdgeBetweenVerticesNC(complex,v1,v2);
    end
);


InstallMethod( OtherVertexOfEdgeNC,
    "for a polygonal complex, a vertex, and an edge",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function( complex, vertex, edge )
        local possVert;
        
        possVert := VerticesOfEdges(complex)[edge];
        if vertex = possVert[1] then
            return possVert[2];
        else
            return possVert[1];
        fi;
    end
);
InstallMethod( OtherVertexOfEdgeNC,
    "for a bend polygonal complex, a vertex, and an edge",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt],
    function( complex, vertex, edge )
        local possVert;
        
        possVert := VerticesOfEdges(complex)[edge];
        if Length(possVert) = 1 then
            return possVert[1];
        fi;

        if vertex = possVert[1] then
            return possVert[2];
        else
            return possVert[1];
        fi;
    end
);
InstallMethod( OtherVertexOfEdge,
    "for a VEF-complex, a vertex, and an edge",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function( complex, vertex, edge )
        local name;
        
        name := "OtherVertexOfEdge";
        __SIMPLICIAL_CheckEdge(complex, edge, name);
        __SIMPLICIAL_CheckIncidenceVertexEdge(complex, vertex, edge, name);
        return OtherVertexOfEdgeNC(complex, vertex, edge);
    end
);



#######################################
##
##      Moving within faces
##
InstallMethod( EdgesInFaceByVerticesNC,
    "for a VEF-complex, a face, and a set of two vertices",
    [IsVEFComplex, IsPosInt, IsSet],
    function( complex, face, vertSet )
        return Filtered( EdgesOfFaces(complex)[face], e ->
            VerticesOfEdges(complex)[e] = vertSet);
    end
);
InstallMethod( EdgesInFaceByVerticesNC,
    "for a VEF-complex, a face, and a list of two vertices",
    [IsVEFComplex, IsPosInt, IsList],
    function( complex, face, vertList )
        return EdgesInFaceByVerticesNC(complex, face, Set(vertList));
    end
);
InstallMethod( EdgesInFaceByVertices,
    "for a VEF-complex, a face, and a set of two vertices",
    [IsVEFComplex, IsPosInt, IsSet],
    function( complex, face, vertSet )
        __SIMPLICIAL_CheckFace(complex, face, "EdgesInFaceByVertices");
        return EdgesInFaceByVerticesNC(complex, face, vertSet);
    end
);
InstallMethod( EdgesInFaceByVertices,
    "for a VEF-complex, a face, and a list of two vertices",
    [IsVEFComplex, IsPosInt, IsList],
    function( complex, face, vertList )
        return EdgesInFaceByVertices(complex, face, Set(vertList));
    end
);


InstallMethod( EdgeInFaceByVerticesNC, 
    "for a VEF-complex, a face, and a set of two vertices",
    [IsVEFComplex, IsPosInt, IsSet],
    function( complex, face, vertSet )
        local possEdges;

        possEdges := EdgesInFaceByVerticesNC(complex, face, vertSet);
        if Length(possEdges) = 0  or Length(possEdges) > 1 then
            return fail;
        fi;
        return possEdges[1];
    end
);
InstallMethod( EdgeInFaceByVerticesNC,
    "for a VEF-complex, a face, and a list of two vertices",
    [IsVEFComplex, IsPosInt, IsList],
    function( complex, face, vertList )
        return EdgeInFaceByVerticesNC(complex, face, Set(vertList));
    end
);
InstallMethod( EdgeInFaceByVertices,
    "for a VEF-complex, a face, and a set of two vertices",
    [IsVEFComplex, IsPosInt, IsSet],
    function( complex, face, vertSet )
        __SIMPLICIAL_CheckFace(complex, face, "EdgeInFaceByVertices");
        return EdgeInFaceByVerticesNC(complex, face, vertSet);
    end
);
InstallMethod( EdgeInFaceByVertices,
    "for a VEF-complex, a face, and a list of two vertices",
    [IsVEFComplex, IsPosInt, IsList],
    function( complex, face, vertList )
        return EdgeInFaceByVertices(complex, face, Set(vertList));
    end
);


######


InstallMethod( OtherEdgesOfVertexInFaceNC,
    "for a VEF-complex, a vertex, an edge, and a face",
    [IsVEFComplex, IsPosInt, IsPosInt, IsPosInt],
    function( complex, vertex, edge, face )
        local possEdges, poss, verts, res;

        possEdges := EdgesOfFaces(complex)[face];
        res := [];
        for poss in possEdges do
            if poss <> edge then
                verts := VerticesOfEdges(complex)[poss];
                if vertex in verts then # Should this be done manually for polygonal complexes?
                    Add(res, poss);
                fi;
            fi;
        od;

        return res;
    end
);
InstallMethod( OtherEdgesOfVertexInFace,
    "for a VEF-complex, a vertex, an edge, and a face",
    [IsVEFComplex, IsPosInt, IsPosInt, IsPosInt],
    function( complex, vertex, edge, face )
        local name;

        name := "OtherEdgesOfVertexInFace";
        __SIMPLICIAL_CheckVertex(complex, vertex, name);
        __SIMPLICIAL_CheckEdge(complex, edge, name);
        __SIMPLICIAL_CheckFace(complex, face, name);
        __SIMPLICIAL_CheckIncidenceVertexEdge(complex, vertex, edge, name);
        __SIMPLICIAL_CheckIncidenceEdgeFace(complex, edge, face, name);

        return OtherEdgesOfVertexInFaceNC(complex, vertex, edge, face);
    end
);


InstallMethod( OtherEdgeOfVertexInFaceNC,
    "for a polygonal complex, a vertex, an edge, and a face",
    [IsPolygonalComplex, IsPosInt, IsPosInt, IsPosInt],
    function( complex, vertex, edge, face )
        local possEdges, poss, verts;

        possEdges := EdgesOfFaces(complex)[face];
        for poss in possEdges do
            if poss <> edge then
                verts := VerticesOfEdges(complex)[poss];
                if vertex = verts[1] or vertex = verts[2] then
                    return poss;
                fi;
            fi;
        od;

        Error("OtherEdgeOfVertexInFaceNC: No valid edge found.");
    end
);
InstallMethod( OtherEdgeOfVertexInFaceNC,
    "for a bend polygonal complex, a vertex, an edge, and a face",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt, IsPosInt],
    function( complex, vertex, edge, face )
        local possEdges;

        possEdges := OtherEdgesOfVertexInFaceNC(complex,vertex,edge,face);
        if Length(possEdges) <> 1 then
            return fail;
        else
            return possEdges[1];
        fi;
    end
);
InstallMethod( OtherEdgeOfVertexInFace,
    "for a VEF-complex, a vertex, an edge, and a face",
    [IsVEFComplex, IsPosInt, IsPosInt, IsPosInt],
    function( complex, vertex, edge, face )
        local name;

        name := "OtherEdgeOfVertexInFace";
        __SIMPLICIAL_CheckVertex(complex, vertex, name);
        __SIMPLICIAL_CheckEdge(complex, edge, name);
        __SIMPLICIAL_CheckFace(complex, face, name);
        __SIMPLICIAL_CheckIncidenceVertexEdge(complex, vertex, edge, name);
        __SIMPLICIAL_CheckIncidenceEdgeFace(complex, edge, face, name);

        return OtherEdgeOfVertexInFaceNC(complex, vertex, edge, face);
    end
);


###########

InstallMethod( OppositeVertexOfEdgeInTriangleNC,
    "for a polygonal complex, an edge, and a face",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function(complex, edge, face)
        local vFace, vEdge, v;

        vFace := VerticesOfFaces(complex)[face];
        vEdge := VerticesOfEdges(complex)[edge];
        for v in vFace do
            if not v in vEdge then
                return v;
            fi;
        od;

        Error("OppositeVertexOfEdgeInTriangleNC: Illegal input. Please don't use the NC-version.");
    end
);
InstallMethod( OppositeVertexOfEdgeInTriangle,
    "for a polygonal complex, an edge, and a face",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function(complex, edge, face)
        local name, len;

        name := "OppositeVertexOfEdgeInTriangle";
        __SIMPLICIAL_CheckEdge(complex, edge, name);
        __SIMPLICIAL_CheckFace(complex, face, name);
        __SIMPLICIAL_CheckIncidenceEdgeFace(complex, edge, face, name);
        len := Length( EdgesOfFaces(complex)[face] );
        if len <> 3 then
            Error(Concatenation(name, ": Given face", String(face), 
                    " should be a triangle, but has ", String(len), " edges." ));
        fi;

        return OppositeVertexOfEdgeInTriangleNC(complex, edge, face);
    end
);


InstallMethod( OppositeEdgeOfVertexInTriangleNC,
    "for a polygonal complex, a vertex, and a face",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function(complex, vertex, face)
        local eFace, eVertex, e;

        eFace := EdgesOfFaces(complex)[face];
        eVertex := EdgesOfVertices(complex)[vertex];
        for e in eFace do
            if not e in eVertex then
                return e;
            fi;
        od;

        Error("OppositeEdgexOfVertexInTriangleNC: Illegal input. Please don't use the NC-version.");
    end
);
InstallMethod( OppositeEdgeOfVertexInTriangle,
    "for a polygonal complex, a vertex, and a face",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function(complex, vertex, face)
        local name, len;

        name := "OppositeEdgeOfVertexInTriangle";
        __SIMPLICIAL_CheckVertex(complex, vertex, name);
        __SIMPLICIAL_CheckFace(complex, face, name);
        __SIMPLICIAL_CheckIncidenceVertexFace(complex, vertex, face, name);
        len := Length( EdgesOfFaces(complex)[face] );
        if len <> 3 then
            Error(Concatenation(name, ": Given face", String(face), 
                    " should be a triangle, but has ", String(len), " edges." ));
        fi;

        return OppositeEdgeOfVertexInTriangleNC(complex, vertex, face);
    end
);


#######################################
##
##      Moving between faces
##
InstallMethod( IsFacesAdjacentNC, "for a VEF-complex and two faces",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, f1, f2)
        if f1 < f2 then
            return [f1,f2] in FacesOfEdges(complex);
        else
            return [f2,f1] in FacesOfEdges(complex);
        fi;
    end
);
InstallMethod( IsFacesAdjacent, "for a VEF-complex and two faces",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, f1, f2)
        local name;

        name := "IsFacesAdjacent";
        __SIMPLICIAL_CheckFace(complex, f1, name);
        __SIMPLICIAL_CheckFace(complex, f2, name);

        return IsFacesAdjacentNC(complex, f1, f2);
    end
);

InstallMethod( EdgesBetweenFacesNC, "for a VEF-complex and two faces",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, f1, f2)
        local res, e1, e2, e;

        res := [];
        e1 := EdgesOfFaces(complex)[f1];
        e2 := EdgesOfFaces(complex)[f2];
        for e in e1 do
            if e in e2 then
                Add(res,e);
            fi;
        od;

        return res;
    end
);
InstallMethod( EdgesBetweenFaces, "for a VEF-complex and two faces",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, f1, f2)
        local name;

        name := "EdgesBetweenFaces";
        __SIMPLICIAL_CheckFace(complex, f1, name);
        __SIMPLICIAL_CheckFace(complex, f2, name);

        return EdgesBetweenFacesNC(complex, f1, f2);
    end
);

InstallMethod( EdgeBetweenFacesNC, "for a VEF-complex and two faces",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, f1, f2)
        local edges;

        edges := EdgesBetweenFacesNC(complex, f1, f2);
        if Length(edges) = 1 then
            return edges[1];
        else
            return fail;
        fi;
    end
);
InstallMethod( EdgeBetweenFaces, "for a VEF-complex and two faces",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function(complex, f1, f2)
        local name;

        name := "EdgeBetweenFaces";
        __SIMPLICIAL_CheckFace(complex, f1, name);
        __SIMPLICIAL_CheckFace(complex, f2, name);

        return EdgeBetweenFacesNC(complex,f1,f2);
    end
);

############

InstallMethod( NeighbourFacesByEdgeNC,
    "for a polygonal complex, a face, and an edge",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function(complex,face,edge)
        local possFaces;

        possFaces := FacesOfEdges(complex)[edge];
        if Length(possFaces) = 1 then
            return [];
        fi;

        if Length(possFaces) = 2 then
            if possFaces[1] = face then
                return [possFaces[2]];
            else
                return [possFaces[1]];
            fi;
        fi;

        return Difference(possFaces, [face]);
    end
);
InstallMethod( NeighbourFacesByEdgeNC,
    "for a bend polygonal complex, a face, and an edge",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt],
    function(complex, face, edge)
        local possFlags, neigh, faces;

        # Which flags are incident to this face-edge-pair?
        possFlags := Intersection( LocalFlagsOfFaces(complex)[face], LocalFlagsOfEdges(complex)[edge] );

        neigh := Filtered( LocalFlagsOfHalfEdges(complex), p -> possFlags[1] in p )[1];
        neigh := Difference( neigh, [possFlags[1]] );

        faces := List(neigh, n -> FacesOfLocalFlags(complex)[n]);
        return Set(faces);
    end
);

InstallMethod( NeighbourFaceByEdgeNC,
    "for a polygonal complex, a face, and an edge",
    [IsPolygonalComplex, IsPosInt, IsPosInt],
    function( complex, face, edge )
        local possFaces;
        
        possFaces := FacesOfEdges(complex)[edge];
        if Length(possFaces) <> 2 then
            return fail;
        fi;

        if possFaces[1] = face then
            return possFaces[2];
        else
            return possFaces[1];
        fi;
    end
);
InstallMethod( NeighbourFaceByEdgeNC,
    "for a bend polygonal complex, a face, and an edge",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt],
    function(complex, face, edge)
        local poss;

        poss := NeighbourFacesByEdgeNC(complex,face,edge);
        if Length(poss) = 1 then
            return poss[1];
        else
            return fail;
        fi;
    end
);


InstallMethod( NeighbourFaceByEdge,
    "for a VEF-complex, a face, and an edge",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function( complex, face, edge )
        local name;
        
        name := "NeighbourFaceByEdge";
        __SIMPLICIAL_CheckEdge(complex,edge, name);
        __SIMPLICIAL_CheckFace(complex, face, name);
        __SIMPLICIAL_CheckIncidenceEdgeFace(complex, edge, face, name);
        
        return NeighbourFaceByEdgeNC(complex, face, edge);
    end
);
InstallMethod( NeighbourFacesByEdge,
    "for a VEF-complex, a face, and an edge",
    [IsVEFComplex, IsPosInt, IsPosInt],
    function( complex, face, edge )
        local name;
        
        name := "NeighbourFacesByEdge";
        __SIMPLICIAL_CheckEdge(complex,edge, name);
        __SIMPLICIAL_CheckFace(complex, face, name);
        __SIMPLICIAL_CheckIncidenceEdgeFace(complex, edge, face, name);
        
        return NeighbourFacesByEdgeNC(complex, face, edge);
    end
);


