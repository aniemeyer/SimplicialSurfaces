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

#! @Chapter Access to bend polygonal complexes
#! @ChapterLabel AccessBendPolygonalComplex
#! 

####
# No matter what, a chapter should start with a short description of its
# contents, together with necessary references to previous material (if
# someone does not read the chapters in order).
####

#! In section <Ref Sect="PolygonalStructures_bend"/> we introduced 
#! the concept of <K>BendPolygonalComplex</K> that went further than
#! a pure incidence geometry. This chapter describes this additional
#! structure and how to access it.
#! 

#! TODO restructure this chapter after more information about correct access design is obtained

#! @Section Local flags
#! @SectionLabel AccessBend_LocalFlags

#! The main difference in contrast to a polygonal complex is an additional
#! structure on the flags of the bend polygonal complex. Since these flags
#! are different from the flags of the incidence structure, we denote them
#! as <E>local flags</E>, i.e. flags of the single polygons (before any
#! identifications took place).
#!
#! We will identify the set of local flags with a set of positive integers.
#! The attribute <K>LocalFlags</K> returns this set. In addition there are
#! two other structures:
#! <Enum>
#!     <Item>The interaction of the local flags with the incidence geometry,
#!           encoded by the attributes <K>VerticesOfLocalFlags</K> or
#!           <K>LocalFlagsOfVertices</K>.</Item>
#!     <Item>The neighbouring relations between the local flags. Generically
#!           they are given by the attribute <K>LocalFlagVertexPartition</K> 
#!           but if all classes have at most two elements, 
#!           <K>LocalFlagVertexInvolution</K> can be used instead.</Item>
#! </Enum>
#!
 
#! @Description
#! Return the set of local flags of the given bend polygonal complex.
#!
#! @Arguments bendComplex
#! @Returns A set of positive integers
DeclareAttribute( "LocalFlags", [IsBendPolygonalComplex] );


#! @BeginGroup
#! @Description
#! Return lists mapping a local flag (represented by its position
#! in <K>LocalFlags</K>(<A>bendComplex</A>)) to the global
#! vertex/edge/face they belong to.
#!
#! @Returns A list of positive integers
#! @Arguments bendComplex
DeclareAttribute( "VerticesOfLocalFlags", [IsBendPolygonalComplex] );
#! @Arguments bendComplex
DeclareAttribute( "EdgesOfLocalFlags", [IsBendPolygonalComplex] );
#! @Arguments bendComplex
DeclareAttribute( "FacesOfLocalFlags", [IsBendPolygonalComplex] );
#! @EndGroup


#! @BeginGroup
#! @Description
#! Return lists mapping a global vertex/edge/face to the set of 
#! local flags surrounding it.
#!
#! @Returns A list of sets of positive integers
#! @Arguments bendComplex
DeclareAttribute( "LocalFlagsOfVertices", [IsBendPolygonalComplex] );
#! @Arguments bendComplex
DeclareAttribute( "LocalFlagsOfEdges", [IsBendPolygonalComplex] );
#! @Arguments bendComplex
DeclareAttribute( "LocalFlagsOfFaces", [IsBendPolygonalComplex] );
#! @EndGroup


#! @BeginGroup
#! @Description
#! Return the set of partitions of the local flags with regard to the
#! vertex/edge/face-equivalence relation. The local flags are given by
#! their positions in <K>LocalFlags</K>(<A>bendComplex</A>).
#!
#! @Returns A set of sets
#! @Arguments bendComplex
DeclareAttribute("LocalFlagVertexPartition", [IsBendPolygonalComplex]);
#! @Arguments bendComplex
DeclareAttribute("LocalFlagEdgePartition", [IsBendPolygonalComplex]);
#! @Arguments bendComplex
DeclareAttribute("LocalFlagFacePartition", [IsBendPolygonalComplex]);
#! @EndGroup


#! @BeginGroup
#! @Description
#! Return the partitions of the local flags with regard to the
#! vertex/edge/face-equivalence relation as involutions. 
#! The local flags are given by
#! their positions in <K>LocalFlags</K>(<A>bendComplex</A>).
#! 
#! If this is not possible, <K>fail</K> is returned instead.
#!
#! @Returns An involution or <K>fail</K>
#! @Arguments bendComplex
DeclareAttribute("LocalFlagVertexInvolution", [IsBendPolygonalComplex]);
#! @Arguments bendComplex
DeclareAttribute("LocalFlagEdgeInvolution", [IsBendPolygonalComplex]);
#! @Arguments bendComplex
DeclareAttribute("LocalFlagFaceInvolution", [IsBendPolygonalComplex]);
#! @EndGroup


#! @Section Bend faces and edges
#! @SectionLabel AccessBend_BendFacesEdges
#!
#! 

