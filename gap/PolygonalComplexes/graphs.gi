#############################################################################
##
##  SimplicialSurface package
##
##  Copyright 2012-2018
##    Markus Baumeister, RWTH Aachen University
##    Alice Niemeyer, RWTH Aachen University 
##
## Licensed under the GPL 3 or later.
##
#############################################################################

#####
#####
#####   This file contains the implementations for the incidence graphs
#####
#####


#######################################
##
##      Digraphs
##
if IsPackageMarkedForLoading( "Digraphs", ">=0.10.1" ) then
    InstallMethod( IncidenceDigraphsGraph, "for a polygonal complex",
        [IsPolygonalComplex],
        function( complex )
            return Digraph( IncidenceGrapeGraph(complex).graph );
            #TODO is there a better way? This ignores the colouring
        end
    );

    InstallMethod( EdgeDigraphsGraph, "for a polygonal complex",
        [IsPolygonalComplex],
        function(complex)
            local graph;

            # Digraphs can only create graphs with vertices [1..n]
            # Therefore we have to take a subgraph of this graph
            graph := DigraphByEdges( Compacted( VerticesOfEdges(complex) ) );
            return InducedSubdigraph( graph, 
                VerticesAttributeOfComplex(complex) );
        end
    );

    InstallMethod(FaceDigraphsGraph, "for a polygonal surface",[IsPolygonalSurface],
	function(surface)
		local i, diedges, arcs, loops;
		
		arcs := [];
		diedges := Filtered( FacesOfEdges(surface), i->Length(i)=2);
		for i in [1..Length(diedges)] do
			arcs[2*i-1] := diedges[i];
			arcs[2*i] := Reversed(diedges[i]);
		od;
		
		loops:=Filtered( FacesOfEdges(surface), i->Length(i)=1);
		for i in [1..Length(loops)] do
			Add(arcs,[loops[i][1],loops[i][1]]);
		od;
		
		return DigraphByEdges(arcs);
	end
);
fi;
##
##      End Digraphs
##
#######################################




#######################################
##
##      GRAPE
##
if IsPackageMarkedForLoading( "GRAPE", ">=0" ) then
    InstallMethod( IncidenceGrapeGraph, "for a polygonal complex",
        [IsPolygonalComplex],
        function(complex)
 	    local graph, vertices, edges, faces, names, colours, incidence, 
	        trivialAction, maxVert, maxEdge, colourClasses;

            maxVert := VerticesAttributeOfComplex(complex)[NumberOfVertices(complex)];
            maxEdge := Edges(complex)[NumberOfEdges(complex)];
            vertices := ShallowCopy( VerticesAttributeOfComplex(complex) );
            edges := List( Edges(complex), e -> e + maxVert );
            faces := List( Faces(complex), f -> f + maxVert + maxEdge );

            names := Union( vertices, edges, faces);
	    colours := [vertices,edges, faces]; # This gives the actual names
	    incidence := function(x,y)
                if x in vertices and y in edges then
                    return x in VerticesOfEdges(complex)[y-maxVert];
                elif x in edges and y in vertices then
                    return y in VerticesOfEdges(complex)[x-maxVert];
                elif x in edges and y in faces then
                    return x-maxVert in EdgesOfFaces(complex)[y-maxVert-maxEdge];
                elif x in faces and y in edges then
                    return y-maxVert in EdgesOfFaces(complex)[x-maxVert-maxEdge];
                else
		    return false;
	        fi; 
	    end;

	    trivialAction := function( pnt, g )
	        return pnt;
	    end;

	    graph := Graph( Group( () ), names, trivialAction, incidence );
            # Since a grape colouring refers to the graph vertices instead of the names, 
            # we need to redefine the colour classes
            colourClasses := [
                [1 .. Length(vertices)],
                Length(vertices) + [1 .. Length(edges)],
                Length(vertices) + Length(edges) + [1 .. Length(faces)]
            ];

	    return rec( graph := graph, colourClasses := colourClasses );   
        end
    );
    InstallMethod( EdgeGrapeGraph, "for a polygonal complex",
        [IsPolygonalComplex],
        function(complex)
    	    local graph, vertices, names, incidence, trivialAction;

	    vertices := VerticesAttributeOfComplex(complex);

            names := vertices;
	    incidence := function(x,y)
		return Set([x,y]) in VerticesOfEdges(complex);
	    end;

	    trivialAction := function( pnt, g )
		return pnt;
	    end;

	    graph := Graph( Group( () ), names, trivialAction, incidence );

	    return graph;
        end
    );
fi;
##
##      End GRAPE
##
#######################################



#######################################
##
##      NautyTracesInterface
##
if IsPackageMarkedForLoading("NautyTracesInterface", ">=0") then
    InstallMethod( IncidenceNautyGraph, "for a polygonal complex",
        [IsPolygonalComplex],
        function(complex)
            local maxVertex, maxEdge, maxFace, edgeList, colourList, v, e, f,
                colSet, vertexList, verts;

            maxVertex := VerticesAttributeOfComplex(complex)[NumberOfVertices(complex)];
            maxEdge := Edges(complex)[NumberOfEdges(complex)];

            vertexList := ShallowCopy( VerticesAttributeOfComplex(complex) );
            edgeList := [];
            colourList := ListWithIdenticalEntries( NumberOfVertices(complex), 0 );

            for e in Edges(complex) do
                # There are two vertices for each edge
                verts := VerticesOfEdges(complex)[e];
                Append( edgeList, [ [verts[1], maxVertex+e], [verts[2], maxVertex+e] ] );
            od;
            Append(vertexList, Edges(complex) + maxVertex);
            Append(colourList, ListWithIdenticalEntries( NumberOfEdges(complex), 1 ));

            for f in Faces(complex) do
                Add(colourList, 2); # done manually since NumberOfFaces is not necessarily computed
                for e in EdgesOfFaces(complex)[f] do
                    Add( edgeList, [maxVertex + e, maxVertex + maxEdge + f] );
                od;
            od;
            Append(vertexList, Faces(complex) + maxVertex + maxEdge);

            return NautyColoredGraphWithNodeLabels( edgeList, colourList, vertexList );
        end
    );

    InstallMethod( ChamberAdjacencyGraph, "for a twisted polygonal complex",
        [IsTwistedPolygonalComplex],
        function(complex)
            local nrNodes, i, chambersOfLabels, c, edges, cl, twoEdges;

            nrNodes := NumberOfChambers(complex);
            chambersOfLabels := [];
            for i in [1 .. nrNodes] do
                c := Chambers(complex)[i];
                chambersOfLabels[c] := i;
            od;

            # construct the edges
            edges := [];
            Add(edges, List(ZeroAdjacencyClasses(complex), cl -> List(cl, c -> chambersOfLabels[c])));
            Add(edges, List(OneAdjacencyClasses(complex), cl -> List(cl, c -> chambersOfLabels[c])));
            twoEdges := [];
            for cl in TwoAdjacencyClasses(complex) do
                if Length(cl) > 1 then
                    Add(twoEdges, List(cl, c -> chambersOfLabels[c]));
                fi;
            od;
            Add(edges, twoEdges);

            return NautyEdgeColoredGraph(edges, nrNodes);
        end
    );

    InstallMethod( EdgeNautyGraph, "for a polygonal complex",
        [IsPolygonalComplex],
        function(complex)
            return NautyGraphWithNodeLabels( VerticesOfEdges(complex), 
                VerticesAttributeOfComplex(complex) );
        end
    );

    InstallMethod(FaceNautyGraph, "for a polygonal surface",[IsPolygonalSurface],
	function(surface)
		local i, diedges, loops;
			
		diedges:=Filtered(FacesOfEdges(surface),i->Length(i)=2);	
		loops:=Filtered(FacesOfEdges(surface),i->Length(i)=1);
		loops:=List(loops,i->[i[1],i[1]]);
			
		return NautyGraphWithNodeLabels( Concatenation(diedges,loops), 
                Faces(surface) );
	end
    );
fi;
##
##      End NautyTracesInterface
##
#######################################




#######################################
##
##      Isomorphism test and automorphism group
##


## The order of installation describes which of these functions is
## preferred - the last one has highest priority.
InstallMethod( IsIsomorphic, "for two twisted polygonal complexes",
    [IsTwistedPolygonalComplex, IsTwistedPolygonalComplex],
    function(complex1, complex2)
        Error("IsIsomorphic: The package NautyTracesInterface has to be available.");
    end
);
InstallOtherMethod( IsIsomorphic, "for two polygonal complexes",
    [IsPolygonalComplex, IsPolygonalComplex],
    function(complex1, complex2)
        Error("IsIsomorphic for polygonal complexes: One of the packages NautyTracesInterface or GRAPE have to be available.");
    end
);
if SIMPLICIAL_ENABLE_SURFACE_REDISPATCH then
    RedispatchOnCondition( IsIsomorphic, true, 
        [IsTwistedPolygonalComplex, IsTwistedPolygonalComplex],
        [IsPolygonalComplex, IsPolygonalComplex], 0);
fi;

InstallMethod( AutomorphismGroup, "for a twisted polygonal complex", 
    [IsTwistedPolygonalComplex],
    function(complex)
        Error("AutomorphismGroup: The package NautyTracesInterface has to be available.");
    end
);
if IsPackageMarkedForLoading( "GRAPE", ">=0" ) then
    InstallOtherMethod( IsIsomorphic,
        "for two polygonal complexes",
        [IsPolygonalComplex, IsPolygonalComplex],
        function(complex1, complex2)
            local inc1, inc2;

            inc1 := IncidenceGrapeGraph(complex1);
            inc2 := IncidenceGrapeGraph(complex2);
            # We copy the structure fully, so that all components stay mutable
            # (this is necessary for GRAPE to function)
            return IsIsomorphicGraph(
                rec( graph := ShallowCopy(inc1.graph), colourClasses := ShallowCopy(inc1.colourClasses) ),
                rec( graph := ShallowCopy(inc2.graph), colourClasses := ShallowCopy(inc2.colourClasses) )
            );
        end
    );
fi;

if IsPackageMarkedForLoading("NautyTracesInterface", ">=0") then
    InstallMethod( IsIsomorphic, 
        "for two twisted polygonal complexes", 
        [IsTwistedPolygonalComplex, IsTwistedPolygonalComplex],
        function(complex1, complex2)
            return IsomorphismGraphs( 
                ChamberAdjacencyGraph(complex1),
                ChamberAdjacencyGraph(complex2)) <> fail;
        end
    );

    InstallMethod( AutomorphismGroup, "for a twisted polygonal complex", 
        [IsTwistedPolygonalComplex],
        function(complex)
            return AutomorphismGroup( ChamberAdjacencyGraph(complex) );
        end
    );
fi;


InstallMethod( IsomorphismRepresentatives,
    "for a list of twisted polygonal complexes", [IsList],
    function(ls)
        local newList, p, q, newOne;

        for p in ls do
            if not IsTwistedPolygonalComplex(p) then
                Error("IsomorphismRepresentatives: Argument has to be a list of twisted polygonal complexes.");
            fi;
        od;

        newList := [];
        for p in ls do
            newOne := true;
            for q in newList do
                if IsIsomorphic(p,q) then
                    newOne := false;
                    break;
                fi;
            od;
            if newOne then
                Add(newList, p);
            fi;
        od;

        return newList;
    end
);

InstallMethod( CanonicalRepresentativeOfPolygonalSurface, 
    "for a polygonal surface", [IsPolygonalSurface],
    function( surf )
        local originalfacesofsurf, originaledgesofsurf, originalverticesofsurf,
        totalgraphverts, mapfaces, mapedges, mapvertices, currentvert, i, vertsofgraph,
        edges, edgesofface, j, verticesofface, verticesofedge, graph, perm, perminv, 
        edgesoffacesofsurf, F, edgesofface2, verticesofedgesofsurf, e, verticesofedge2,
        newfaces, newedges, newvertices, n1, n2, n3,
        mapfaces2, mapedges2, mapvertices2, edgesoffacesofsurf2, verticesofedgesofsurf2,
        surf2, surf3, colours, inversefacemap, inverseedgemap, inversevertexmap;


        # The original faces/edges/vertices of surf
        originalfacesofsurf := Faces(surf); 
        originaledgesofsurf := Edges(surf);
        originalverticesofsurf := VerticesAttributeOfComplex(surf);
        # The number of faces/edges/vertices of surf
        n1 := NumberOfFaces(surf);
        n2 := NumberOfEdges(surf);
        n3 := NumberOfVertices(surf);
        # Total number of elements of surf - we make a bijection with the elements
        totalgraphverts := n1+n2+n3;

        # Create maps from the elements to the new labels.
        # Map faces to [1 .. n1], edges to [n1+1 .. n1+n2] and vertices to [n1+n2+1 .. totalverts]
        # Let i be a face, then the i maps to mapfaces[i]
        mapfaces := [];
        mapedges := [];
        mapvertices := [];
        # Also assign each element a colour - faces are 1, edges are 2, vertices are 3
        # Let i be an element of surf with the new labelling. Then the colour can be
        # established with colours[i]
        colours:=[];
        for i in [1..n1] do
            mapfaces[originalfacesofsurf[i]] := i;
        od;
        Append(colours, ListWithIdenticalEntries(n1,1));
        for i in [1..n2] do
            mapedges[originaledgesofsurf[i]] := i + n1;
        od;
        Append(colours, ListWithIdenticalEntries(n2,2));
        for i in [1..n3] do
            mapvertices[originalverticesofsurf[i]] := i + n1 + n2;
        od;
        Append(colours, ListWithIdenticalEntries(n3,3));


        # Create the corresponding graph.
        # The elements of surf form the vertices of the graph.
        # There is an edge in the graph if the corresponding elements of surf are incident.
        vertsofgraph := [1 .. totalgraphverts];

        # Loop through the faces and for each face loop through vertices and edges.
        # If the face is incident with the edge/vertex, put an edge in the graph.
        edges := [];
        for i in originalfacesofsurf do
            edgesofface := EdgesOfFaces( surf )[i];
            for j in edgesofface do
                Add(edges, [mapfaces[i], mapedges[j]]);
            od; 
            verticesofface := VerticesOfFaces( surf ) [i]; 
            for j in verticesofface do
                Add(edges, [mapfaces[i], mapvertices[j]]);
            od; 
        od;
        # Repeat for edges. Only need to check for vertices, since we did faces-edges above.
        for i in originaledgesofsurf do
            verticesofedge := VerticesOfEdges( surf )[i];
            for j in verticesofedge do
                Add(edges, [mapedges[i], mapvertices[j]]);
            od; 
        od;
        # No need to loop through vertices, since vertices are not incident with vertices,
        # and incidence with faces/edges done above.

        # Find the canonical form of the graph with Nauty,
        # preserving the colours (ie  map faces to faces, edges to edge etc...)
        graph := NautyColoredGraph(edges, colours);
        # Find the permutation which will map the old vertices to the new.
        perminv := CanonicalLabeling(graph);
        perm := perminv^(-1);

        # Now that we have the canonical labelling, we can write the surface in canonical form.
        # First take a face and write it by its edges
        # Map the face to its new labelling, and then permute to canonical form.
        # Map the edges to the new labelling
        # Permute the newly labelled eges to their canonical form.
        # Now put the edges in the list of FacesByEdges
        edgesoffacesofsurf := [];
        for i in originalfacesofsurf do
            F := mapfaces[i]^perm;
            edgesofface := EdgesOfFaces( surf )[i];
            edgesofface2 := [];
            for j in [1..Length(edgesofface)] do
                edgesofface2[j] := mapedges[edgesofface[j]]^perm;
            od;
            edgesoffacesofsurf[F] := edgesofface2;
        od;

        # Same as above, but for edges with respect to vertices
        verticesofedgesofsurf := [];
        for i in originaledgesofsurf do
            e := mapedges[i]^perm;
            verticesofedge := VerticesOfEdges( surf )[i];
            verticesofedge2 := [];
            for j in [1..Length(verticesofedge)] do
                verticesofedge2[j] := mapvertices[verticesofedge[j]]^perm;
            od;
            verticesofedgesofsurf[e] := verticesofedge2;
        od;


        # Map the elements to the new labelling, then permute to canonical form
        newfaces := [1..n1];
        newedges := [n1+1..n1+n2];
        newvertices := [n1+n2+1..n1+n2+n3];

        # Now that we have the newly labelled, canonical elements,
        # we map them to a lex least labelling.
        # The newfaces are [1..n1] and will stay that way
        # The newedges are [n1+1..n1+n2] and their minimal labeling is [1..n2]
        #    (so the performed operation is subtraction of n1)
        # newvertices are [n1+n2+1..n1+n2+n3] and their minimal labeling is [1..n3]
        #    (so the performed operation is subtraction of n1+n2)

        edgesoffacesofsurf2 := [];
        for i in newfaces do
            # we talk about the face i
            edgesoffacesofsurf2[i] := [];
            for j in edgesoffacesofsurf[i] do
                Add(edgesoffacesofsurf2[i], j-n1); # shift edge labels by n1
            od;
        od;

        verticesofedgesofsurf2 := [];
        for i in newedges do
            # we talk about the edge i-n1
            verticesofedgesofsurf2[i-n1] := [];
            for j in verticesofedgesofsurf[i] do
                Add(verticesofedgesofsurf2[i-n1], j - n1 - n2); # shift vertex labels by n1+n2
            od;
        od;

        # To get the inverse map (from the new labels to the old ones), 
        # we reverse the lex least map,
        # then the canonical permuting, then the bijection to new labelling.
        inversefacemap := [];
        for i in newfaces do
            inversefacemap[i] := Faces(surf)[i^perminv];
        od;
        inverseedgemap := [];
        for i in newedges do
            inverseedgemap[i-n1] := Edges(surf)[i^perminv - n1];
        od;
        inversevertexmap := [];
        for i in newvertices do
            inversevertexmap[i-n1-n2] := VerticesAttributeOfComplex(surf)[i^perminv - n1 - n2];
        od;

        surf2 := PolygonalSurfaceByDownwardIncidenceNC(verticesofedgesofsurf2, edgesoffacesofsurf2);

        # return the canonical form of the surface and
        # the bijections mapping the new elements to old, by element i in canonical surface
        # maps to inversemap[i] in the old surface.
        return [surf2, PolygonalMorphismByListsNC(surf2, surf, 
                    inversevertexmap, inverseedgemap, inversefacemap)];
    end
);
if SIMPLICIAL_ENABLE_SURFACE_REDISPATCH then
    RedispatchOnCondition( CanonicalRepresentativeOfPolygonalSurface, true, 
        [IsPolygonalComplex], [IsPolygonalSurface], 0 );
fi;


#######################################
##
##      Automorphism group
##

BindGlobal( "__SIMPLICIAL_RestrictToVertices",
    function(complex, g)
        local maxVert, permList, c, vOfC;

        maxVert := Maximum(Vertices(complex));
        permList := [1..maxVert];
        vOfC := VerticesOfChambers(complex);
        for c in Chambers(complex) do
            permList[vOfC[c]] := vOfC[c^g];
        od;
        return PermList(permList);
    end
);
BindGlobal( "__SIMPLICIAL_RestrictToEdges",
    function(complex,  g)
        local maxVert, permList, c, eOfC;

        maxVert := Maximum(Vertices(complex));
        permList := [1..maxVert];
        eOfC := EdgesOfChambers(complex);
        for c in Chambers(complex) do
            permList[eOfC[c]] := eOfC[c^g];
        od;
        return PermList(permList);
    end
);

BindGlobal( "__SIMPLICIAL_RestrictToFaces",
    function(complex, g)
        local maxVert, permList, c, fOfC;

        maxVert := Maximum(Vertices(complex));
        permList := [1..maxVert];
        fOfC := FacesOfChambers(complex);
        for c in Chambers(complex) do
            permList[fOfC[c]] := fOfC[c^g];
        od;
        return PermList(permList);
    end
);

InstallMethod( DisplayAsAutomorphism, 
    "for a polygonal complex and a permutation",
    [IsTwistedPolygonalComplex, IsPerm],
    function(complex, perm)
        local autVert, autEdge, autFace;

        autVert := __SIMPLICIAL_RestrictToVertices(complex, perm);
        autEdge := __SIMPLICIAL_RestrictToEdges(complex, perm );
        autFace := __SIMPLICIAL_RestrictToFaces(complex, perm);

        return [autVert, autEdge, autFace];
    end
);

InstallMethod( AutomorphismGroupOnVertices, "for a polygonal complex",
    [IsTwistedPolygonalComplex],
    function(complex)
        local aut, grp, gens, g;

        aut := AutomorphismGroup(complex);
        gens := [];
        for g in GeneratorsOfGroup(aut) do
            Add( gens, __SIMPLICIAL_RestrictToVertices(complex, g) );
        od;

        return Group( gens );
    end
);
InstallMethod( AutomorphismGroupOnEdges, "for a polygonal complex",
    [IsTwistedPolygonalComplex],
    function(complex)
        local aut, grp, gens, g;

        aut := AutomorphismGroup(complex);
        gens := [];
        for g in GeneratorsOfGroup(aut) do
            Add( gens, __SIMPLICIAL_RestrictToEdges(complex, g) );
        od;

        return Group( gens );
    end
);
InstallMethod( AutomorphismGroupOnFaces, "for a polygonal complex",
    [IsTwistedPolygonalComplex],
    function(complex)
        local aut, grp, gens, g;

        aut := AutomorphismGroup(complex);
        gens := [];
        for g in GeneratorsOfGroup(aut) do
            Add( gens, __SIMPLICIAL_RestrictToFaces(complex, g) );
        od;

        return Group( gens );
    end
);


InstallMethod( IsAutomorphismDefinedByVertices, "for a polygonal complex",
    [IsTwistedPolygonalComplex],
    function(complex)
        return Length( Set(VerticesOfEdges(complex)) ) = NumberOfEdges(complex) 
            and Length( Set(VerticesOfFaces(complex)) ) = NumberOfFaces(complex);
    end
);
InstallMethod( IsAutomorphismDefinedByEdges, "for a polygonal complex",
    [IsTwistedPolygonalComplex],
    function(complex)
        return Length( Set(EdgesOfVertices(complex)) ) = NumberOfVertices(complex) 
            and Length( Set(EdgesOfFaces(complex)) ) = NumberOfFaces(complex);
    end
);
InstallMethod( IsAutomorphismDefinedByFaces, "for a polygonal complex",
    [IsTwistedPolygonalComplex],
    function(complex)
        return Length( Set(FacesOfEdges(complex)) ) = NumberOfEdges(complex) 
            and Length( Set(FacesOfVertices(complex)) ) = NumberOfVertices(complex);
    end
);


InstallMethod( OnVertexEdgePaths, 
    "for a vertex-edge-path on a polygonal complex and an automorphism",
    [IsVertexEdgePath , IsPerm],
    function( vePath, aut )
        local surf, list, display, verts, edges, i;

        surf := AssociatedPolygonalComplex(vePath);
        display := DisplayAsAutomorphism(surf, aut);
        if display = fail then
            return fail;
        fi;

        list := [];
        verts := VerticesAsList(vePath);
        edges := EdgesAsList(vePath);
        for i in [1..Length(verts)] do
            list[2*i-1] := verts[i]^display[1];
        od;
        for i in [1..Length(edges)] do
            list[2*i] := edges[i]^display[2];
        od;

        return VertexEdgePathNC(surf, list);
    end
);


InstallMethod( OnEdgeFacePaths, 
    "for an edge-face-path on a polygonal complex and an automorphism",
    [IsEdgeFacePath , IsPerm],
    function( efPath, aut )
        local surf, list, display, faces, edges, i;

        surf := AssociatedPolygonalComplex(efPath);
        display := DisplayAsAutomorphism(surf, aut);
        if display = fail then
            return fail;
        fi;

        list := [];
        edges := EdgesAsList(efPath);
        faces := FacesAsList(efPath);
        for i in [1..Length(edges)] do
            list[2*i-1] := edges[i]^display[2];
        od;
        for i in [1..Length(faces)] do
            list[2*i] := faces[i]^display[3];
        od;

        return EdgeFacePathNC(surf, list);
    end
);

