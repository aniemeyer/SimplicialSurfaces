#############################################################################
##
##  SimplicialSurface package
##
##  Copyright 2012-2016
##    Markus Baumeister, RWTH Aachen University
##    Alice Niemeyer, RWTH Aachen University 
##
## Licensed under the GPL 3 or later.
##
#############################################################################


BindGlobal( "__SIMPLICIAL_ComputeLocalFlags", 
    function( faces, sizesOfFaces )
        local res, f, i, max;

        res := [];
        for f in faces do
            max := sizesOfFaces[f]; 
            for i in [1..max] do
                if i = 1 then
                    Add(res, [i,max]);
                else
                    Add(res, [i,i-1]);
                fi;

                if i = max then
                    Add(res, [i,1]);
                else
                    Add(res, [i,i+1]);
                fi;
            od;
        od;

        return res;
    end
);


BindGlobal("__SIMPLICIAL_InvolutionFromPartition",
    function(partition)
        local inv, pair;

        inv := ();
        for pair in partition do
            if Length(pair) = 2 then
                inv := inv * (pair[1],pair[2]);
            elif Length(pair) > 2 then
                return fail;
            fi;
        od;
        return inv;
    end
);

InstallMethod( LocalFlagVertexInvolution, "for a bend polygonal complex",
    [IsBendPolygonalComplex],
    function(bendComplex)
        return __SIMPLICIAL_InvolutionFromPartition(
            LocalFlagsOfLocalVertices(bendComplex) );
    end
);
InstallMethod( LocalFlagEdgeInvolution, "for a bend polygonal complex",
    [IsBendPolygonalComplex],
    function(bendComplex)
        return __SIMPLICIAL_InvolutionFromPartition(
            LocalFlagsOfLocalEdges(bendComplex) );
    end
);
InstallMethod( LocalFlagFaceInvolution, "for a bend polygonal complex",
    [IsBendPolygonalComplex],
    function(bendComplex)
        return __SIMPLICIAL_InvolutionFromPartition(
            LocalFlagsOfHalfEdges(bendComplex) );
    end
);


#TODO give them checks and NC-versions
InstallMethod( IsLocalFlagsVertexEquivalent,
    "for a bend polygonal complex and two local flags",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt],
    function(bendComplex, flag1, flag2)
        return flag1^LocalFlagVertexInvolution(bendComplex) = flag2;
    end
);
InstallMethod( IsLocalFlagsEdgeEquivalent,
    "for a bend polygonal complex and two local flags",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt],
    function(bendComplex, flag1, flag2)
        return flag1^LocalFlagEdgeInvolution(bendComplex) = flag2;
    end
);
InstallMethod( IsLocalFlagsFaceEquivalent,
    "for a bend polygonal complex and two local flags",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt],
    function(bendComplex, flag1, flag2)
        local inv, part, bool1, bool2;

        inv := LocalFlagFaceInvolution(bendComplex);
        if inv = fail then
            for part in LocalFlagsOfHalfEdges(bendComplex) do
                bool1 := flag1 in part;
                bool2 := flag2 in part;
                if bool1 and bool2 then
                    return true;
                fi;
                if bool1 or bool2 then
                    return false;
                fi;
            od;
            return false;
        else
            return flag1^LocalFlagFaceInvolution(bendComplex) = flag2;
        fi;
    end
);

InstallMethod( LocalFlagByLocalVertexLocalEdgeFace,
    "for a bend polygonal complex, a local vertex, a local edge and a face",
    [IsBendPolygonalComplex, IsPosInt, IsPosInt, IsPosInt],
    function(bendComplex, localVertex, localEdge, face)
        local flags, locVerts, locEdges, faces, i;

        flags := LocalFlags(bendComplex);
        locVerts := LocalVerticesOfLocalFlags(bendComplex);
        locEdges := LocalEdgesOfLocalFlags(bendComplex);
        faces := FacesOfLocalFlags(bendComplex);
        for i in flags do
            if locVerts[i] = localVertex and locEdges[i] = localEdge and faces[i] = face then
                return i;
            fi;
        od;
        return fail;
    end
);

InstallMethod( BendPolygonalComplexBySignedFacePerimeters,
    "for a list of signed face perimeters", [IsList],
    function(signedFacePerimeter)
        local currLength, vertexMap, edgeMap, edgeSignMap, faceMap,
            vertexPart, edgePart, facePart, edgeOfFace, vertexOfEdge,
            i, size, perimeter, edges, f, edge, verts, found, p, complex;

        currLength := 0;

        vertexMap := []; # VerticesOfLocalFlags
        edgeMap := [];  # EdgesOfLocalFlags
        edgeSignMap := [];
        faceMap := [];  # FacesOfLocalFlags

        vertexPart := []; # LocalFlagsOfLocalVertices
        edgePart := []; # LocalFlagsOfLocalEdges
        facePart := []; # LocalFlagsOfHalfEdges

        edgeOfFace := [];
        vertexOfEdge := [];
        for i in [1..Length(signedFacePerimeter)] do
            if IsBound(signedFacePerimeter[i]) then
                size := Length(signedFacePerimeter[i]);
                # new flags: currLength+1,...,currLength+size
                perimeter := signedFacePerimeter[i];
                Add(perimeter, perimeter[1]);

                edges := [];

                for f in [1..size] do
                    if IsEvenInt(f) then
                        vertexMap[currLength+f] := perimeter[f+1];
                        edgeMap[currLength+f] := AbsoluteValue(perimeter[f]);
                        edgeSignMap[currLength+f] := SignInt(perimeter[f]);
                        if f = size then
                            Add(vertexPart, [currLength+1, currLength+size]);
                        else
                            Add(vertexPart, [currLength+f, currLength+f+1]);
                        fi;
                    else
                        vertexMap[currLength+f] := perimeter[f];
                        edge := AbsoluteValue(perimeter[f+1]);
                        edgeMap[currLength+f] := edge;
                        edgeSignMap[currLength+f] := SignInt(perimeter[f+1])*(-1);
                        Add(edgePart, [currLength+f, currLength+f+1]);

                        Add(edges, edge);
                        verts := Set([perimeter[f],perimeter[f+2]]);
                        if IsBound(vertexOfEdge[edge]) then
                            if vertexOfEdge[edge] <> verts then
                                Error("Inconsistent data");
                            fi;
                        else
                            vertexOfEdge[edge] := verts;
                        fi;
                    fi;
                    faceMap[currLength+f] := i;
                od;

                edgeOfFace[i] := Set(edges);
                currLength := currLength + size;
            fi;
        od;

        # Generate the face partition
        for i in [1..currLength] do
            found := false;
            for p in facePart do
                if edgeMap[i] = edgeMap[p[1]] and edgeSignMap[i] = edgeSignMap[p[1]] then
                    found := true;
                    Add(p,i);
                fi;
            od;
            if not found then
                Add(facePart, [i]);
            fi;
        od;


        complex := Objectify(BendPolygonalComplexType, rec());
        SetVerticesOfEdges(complex, vertexOfEdge);
        SetEdgesOfFaces(complex, edgeOfFace);
        
        SetVerticesOfLocalFlags(complex, vertexMap);
        SetEdgesOfLocalFlags(complex, edgeMap);
        SetFacesOfLocalFlags(complex, faceMap);
        
        SetLocalFlagsOfLocalVertices(complex, vertexPart);
        SetLocalFlagsOfLocalEdges(complex, edgePart);
        SetLocalFlagsOfHalfEdges(complex, facePart);
        return complex;
    end
);
