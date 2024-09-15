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

# This file contains tests for methods to generate surface via Butterfly insertion

BindGlobal( "__SIMPLICIAL_Test_SurfaceGeneration", function()

        local surf, newsurf, t, disc, allsurfs;

        # test the One Face
        surf := SimplicialSurfaceByVerticesInFaces( [[2,3,5] ]);
        umbdescriptor := UmbrellaDescriptorOfSurface(surf);
        newsurf :=  SimplicialSurfaceByUmbrellaDescriptor(umbdescriptor);
        SIMPLICIAL_TestAssert( IsIsomorphic(surf,newsurf) );


       disc := SimplicialSurfaceByDownwardIncidence(
[ , [ 27, 36 ],, [ 44, 57 ],,,, [ 44, 56 ],,,,, [ 35, 47 ],, [ 35, 46 ], [ 36, 46 ],,,,,,,,,, [ 27, 59 ], [ 27, 56 ], [ 27, 46 ], [ 46, 59 ],,,, 
  [ 33, 56 ], [ 33, 35 ], [ 33, 36 ], [ 35, 36 ], [ 36, 56 ],,,,,, [ 44, 47 ], [ 47, 58 ], [ 44, 58 ], [ 46, 58 ], [ 46, 47 ],,,,,,,,, [ 57, 59 ], [ 56, 57 ], [ 57, 58 ], 
  [ 58, 59 ], [ 56, 59 ] ],
[ [ 2, 16, 28 ], [ 2, 27, 37 ], [ 4, 8, 57 ], [ 4, 45, 58 ], [ 29, 46, 59 ],,,,,, [ 15, 16, 36 ], [ 13, 15, 47 ],,,,,,,,,, [ 26, 27, 60 ], [ 26, 28, 29 ],,,,,,,, 
  [ 34, 35, 36 ], [ 33, 35, 37 ],,,,,,,, [ 43, 44, 45 ], [ 44, 46, 47 ],,,,,,,,,,, [ 56, 57, 60 ], [ 56, 58, 59 ] ]);

        # this path goes via a boundary vertex
        t := [ 36, 56, 57 ];
        newsurf := ButterflyInsertion(disc,t)[1];
        SIMPLICIAL_TestAssert( ListCounter(CounterOfVertices(newsurf))=
           [[2,1],[3,3],[4,2],[5,4],[6,2]] and IsEssentialDisc(newsurf));

        
        # test the double Janus head, which has 2 ears.
        surf := SimplicialSurfaceByDownwardIncidence( 
                    [[1,3],[2,3],[1,2],[1,2],[1,4],[2,4]],
                    [[1,2,3],[1,2,4],[3,5,6],[4,5,6]] );
         # this path has inner edges and vertices
         t := [3,1,4];
         newsurf := ButterflyInsertion(surf,t)[1];
         SIMPLICIAL_TestAssert( ListCounter(CounterOfVertices(newsurf))=
        [ [ 3, 2 ], [ 4, 3 ] ] and EulerCharacteristic(newsurf)=2 );

         allsurfs:=AllSimplicialSurfacesByEssentialButterflyInsertion(newsurf);
         SIMPLICIAL_TestAssert( Length(allsurfs)=2 and
        Length(Faces(allsurfs[1])) = 8 and Length(Faces(allsurfs[2])) = 8 );

         
end;
