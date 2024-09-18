

#############################################################################
##
#F IsEssentialDisc( <disc> ) . . . . . . . . test whether <disc> is essential
##
## A simplicial surface <disc> is called an essential simplicial disc, if
## it is a connected simplicial surface of Euler Characteristic 1, no
## inner vertex has degree less than 4, no boundary vertex has degree 1,
## no inner edge connects two boundary vertices and the disc has no 3-waists.
##

InstallMethod( IsEssentialDisc, 
    "for a simplicial surface",
    [IsTwistedPolygonalComplex],
    
      function( disc )

           local bound, v, inner, innerE, e;

           if not IsSimplicialSurface(disc) then return false; fi;

           if not IsConnected(disc) then return false; fi;
           if EulerCharacteristic(disc)<>1 then return false; fi;

           bound := BoundaryVertices(disc);
           # ensure that all boundary vertices have degree at least 2
           for v in bound do
               if DegreeOfVertex(disc,v) < 2 then return false; fi;
           od;
           inner := InnerVertices(disc);
           # ensure that all inner vertices have degree at least 4
           for v in inner do
               if DegreeOfVertex(disc,v) <= 3 then return false; fi;
           od;
           innerE := InnerEdges(disc);
            # ensure that all inner edges do not have two boundary vertices
           for e in innerE do
               v := VerticesOfEdge(disc,e);
               if  v[1] in bound and v[2] in bound then return false; fi;
           od;
           # ensure that the disc has no 3-waists
           if Length(AllThreeWaistsOfComplex(disc)) > 0 then return false; fi;
           return true;
end
);



#############################################################################
##
## A relevant 2-path in an essential disc D is a pair of edges of D sharing
## exactly one vertex W such that the two edges are not adjacent to one
## common face.
## Find all relevant 2-paths up to the action of the automorphism
## group of the surface <surf>
##

BindGlobal( "__SIMPLICIAL_AllEssentialTwoPaths",
    function( surf )

        local vertices, i, t, umbti, u, relpaths, mp, twosets, isgood,
              aut, orbs, isNewInOrb, li;
        umbti := UmbrellaTipDescriptorOfSurface(surf);
        relpaths := [];

        #check if a two-set t contains two neighbouring vertices
        isgood := function(t,u)
            if IsPerm(u) then
                if t[1]^u=t[2] then return false; fi;
                if t[2]^u=t[1] then return false; fi;
            else
               if Position(u,t[1]) = Position(u,t[2]) + 1 or
                  Position(u,t[1]) = Position(u,t[2]) - 1 then
                    return false;
               fi;
            fi;
            return true;
        end;

        isNewInOrb := function(t)

            local j;

            if Length(orbs) = 0 then return true; fi;

            for j in [1 .. Length(orbs)] do
                if t in orbs[j] then return false; fi;
            od;

            return true;
        end;


        aut := AutomorphismGroupOnVertices(surf);
        orbs := [];

        for i in [1 .. Length(umbti)] do
            if IsBound(umbti[i]) then
                # Construct all 2-paths with middle vertex i
                u := umbti[i];
                if IsPerm(u) then
                    # vertex i is inner
                    mp := MovedPoints(u);
                else
                    # vertex i is a boundary vertex. 
                    mp := umbti[i];
                fi;
                # from every vertex in mp we can find a 2-path
                # to any other via i. Just remove non-relevant paths
                twosets := Combinations(mp, 2);
                twosets := Filtered(twosets, t-> isgood(t,u));
                li := List(twosets, t-> [t[1],i,t[2]]);
                for t in li do
                    if isNewInOrb(t) then
                        Add(orbs, Orbit( aut, t, OnTuples) );
                        Add(relpaths, t);
                    fi;
                od;
             fi;
        od;



        return relpaths;

end
);







# Compute up to isomorphism all surfaces that can be obtained from
# the surface surf by a butterfly insertion along a relevant 2-path

InstallMethod( AllSimplicialSurfacesByEssentialButterflyInsertion,
    "for a simplicial surface",
    [IsSimplicialSurface],
    

    function(surf)

        local allp, surfaces;

        allp := __SIMPLICIAL_AllEssentialTwoPaths(surf);

        surfaces := List(allp, t-> ButterflyInsertion(surf, t)[1]);


        return IsomorphismRepresentatives(surfaces);

end);


# All Simplicial discs on <nrFaces> faces and with boundary length <bdLen>.
InstallMethod( AllSimplicialEssentialDiscs,
    "for a pair of positive integers",
    [IsPosInt,IsPosInt],
    
       function( nrFaces, bdLen)

        local reps, bgon, n, newsurfs, surfaces, surf, allp, canrep;

        reps := [];

        canrep := function(sisu)
            local rep;
            rep := CanonicalRepresentativeOfPolygonalSurface(sisu);
            if rep=fail then return fail; fi;
            return rep[1];
        end;


        bgon := canrep(SimplicialUmbrella(bdLen));
        # The simplicial umbrella has as many faces as its
        # boundary length
        n := bdLen;
     
        if nrFaces < n or nrFaces mod 2 <> n mod 2 then
            return [];
        fi;

        reps := [bgon];
        while n < nrFaces do
            newsurfs := Set([]);
            for surf in reps do
                allp := __SIMPLICIAL_AllEssentialTwoPaths(surf);
                surfaces := List(allp, t-> canrep(
                                 ButterflyInsertion(surf, t)[1]));
                newsurfs := Union(newsurfs,surfaces);
            od;
            n := n+2;
            reps := newsurfs;
        od;

        return reps;

end
);




InstallMethod( IsomorphismRepresentativesOfZippedDiscs,

    "for a pair of simplicial discs",
    [IsSimplicialSurface,IsSimplicialSurface],
    
    function(disc1,disc2)

    local  j, perm, spheres, path1, path2, dyclet, n, s, vert1,
           dih, grp1, grp2, double, bound2;


        if IsClosedSurface(disc1) or IsClosedSurface(disc2) or
        not IsConnected(disc1) or not IsConnected(disc2) or
        not EulerCharacteristic(disc1) = 1 or 
        not EulerCharacteristic(disc2) = 1  then
            ErrorNoReturn("both surfaces must be simplicial discs");
        fi;;

        path1 := PerimeterOfHoles(disc1)[1];
        bound2 := PerimeterOfHoles(disc2)[1];
        vert1 := VerticesAsList(path1);
        dyclet := VerticesAsList(bound2);

        if Length(vert1) <> Length(dyclet) then
            ErrorNoReturn("Boundaries must have same length");
        fi;


        grp1 := AutomorphismGroupOnVertices(disc1);
        grp2 := AutomorphismGroupOnVertices(disc2);
        # consider the induced action on the boundary
        grp1 := Action(grp1, vert1{[1..Length(vert1)-1]}, OnPoints);
        grp2 := Action(grp2, dyclet{[1..Length(dyclet)-1]}, OnPoints);

        dyclet := dyclet{[1..Length(dyclet)-1]};
        dih := DihedralGroup( IsPermGroup, 2*Length(dyclet) );

        # compute the double costs in the dihedral group
        # We act with G on U from the right and with H on U from left-inverse
        double :=Orbits(grp1,Elements(dih),OnLeftInverse);
        s:=Orbits(grp2,double,OnRight);
        # double consists of the double cosets grp1\dih/grp2
        double:=Set(List(s,r->Union(r)));

        spheres := [];
        
        n := Length(dyclet);
        # Generate all cyclic permutations and their mirror images
        for j in [1..Length(double)] do
            perm := Permuted(dyclet, double[j][1]);
            perm[Length(perm)+1] := perm[1];
            path2 := VertexEdgePathByVertices(disc2,perm);
            s := JoinVertexEdgePaths(disc1,path1,disc2,path2);
            Add(spheres, s[1]);
        od;

    return IsomorphismRepresentatives(spheres);
end

);


#############################################################################
##
#F Given a subdisc <disc> of a simplicial sphere <sph> returns the comple-
## menting subdisc
##
ComplementingDiscInSphere := function( sphere, disc)

        local fs, fd;

        fs := Faces(sphere);
        fd := Faces(disc);

        if not IsSubset(fs,fd) then
            ErrorNoReturn("disc is not a subdisc of sphere");
        fi;
        return SubsurfaceByFaces(sphere,Difference(fs,fd));

end;


AllDecompositionsSphereAut := function( n )

        local discs, numbers, k, sn_k, y, i, j, s, aut1, aut2;

        if not IsEvenInt(n) then Error("n has to be even"); fi;
        n := n/2;
        # n is the size of a hemisphere

        discs := [];
        numbers := [];
#        for k in [ 4 .. n ] do
        for k in [ 6..7 ] do
            # find zippings of a sphere with <n> faces into
            # two discs whose boundary length is <k>
            sn_k := AllSimplicialEssentialDiscs(n,k);
            y := [];
            for i in [1..Length(sn_k)] do
                aut1:= AutomorphismGroupOnVertices(sn_k[i]);
                for j in [i..Length(sn_k)] do
                    aut2:= AutomorphismGroupOnVertices(sn_k[j]);
                    s := IsomorphismRepresentativesOfZippedDiscs(
                        sn_k[i],sn_k[j]);
                    Print(" k = ", k, " # ", Length(s), "\n");
                    Append(y,s);
                od;
            od;
            y := IsomorphismRepresentatives(y);
            if y <> [] then
                discs[k] := y;
                numbers[k] := Length(y);
            fi;
        od;

        return [discs, numbers];
end;





