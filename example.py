""" example of ofblockmeshdicthelper
try to generate wedged pype object shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""

from __future__ import unicode_literals, print_function, division
import math
from ofblockmeshdicthelper import Vertex, BlockMeshDict, rotate, mirror, \
        translate, SimpleGrading

# Define parameters
wd = wedgedegree = 5.0
radius_x = 0.19/2 # Diameter is 0.19 in the wiki page
length_z = 1.1

# Definition of vertices
v=['',]*8
v[0] = Vertex(0,0,0) # can be written as V[0] = Vertex((0,0,0)).
# Create a new vertex by rotation
v[1] = rotate(Vertex(radius_x, 0, 0), \
        center=v[0], axis=(0,0,1), angle= math.radians(-wd/2))
v[3],v[2]=[rotate(vertex,center=v[0],axis=(0,0,1),angle=math.radians(wd)) \
        for vertex in v[0:2]]
# Create a new vertex by translation
v[4:8]=[translate(vertex,vector=(0,0,length_z)) for vertex in v[0:4]]

# Vertex can be defined as either Vertex class or tuple
arc_p1 = Vertex(v[5].x*1.5, v[5].y*1.5, v[5].z/2)
arc_p2 = (v[6].x*1.5, v[6].y*1.5, v[6].z/2)

spline_p11 = Vertex(v[1].x/3, v[1].y/3, v[1].z-0.02)
spline_p12 = Vertex(v[1].x/3*2, v[1].y/3*2, v[1].z+0.02)

# Create a new vertex by reflection
spline_p21 = mirror(spline_p11, center=v[0], normal=(0,1,0))
spline_p22 = mirror(spline_p12, center=v[0], normal=(0,1,0))

# The BlockMeshDict class can be constructed in serveral ways:
# b=BlockMeshDict(a_vertices_list, cells=???, region=???, ...)
# b=BlockMeshDict(a_block, cells=modified_cells,region=modified_region, ...)
# b=BlockMeshDict(another_BlockMeshDict, metric=modified_metric, 
#                 mergePatchPairs = modified_mergePatchPairs)
b1 = BlockMeshDict(v, cells=(5, 1, 75), region="fluid",\
        boundaries={'xm':'empty_axis' , 'xp':'wall_tankWall', 
                    'ym':'wedge_front', 'yp':'wedge_back', \
                    'zm':'patch_inlet', 'zp':'patch_outlet'}, \
        grading=SimpleGrading(0.1,\
                            ((0.2, 0.3, 4),\
                            (0.6, 0.4, 1),\
                            (0.2, 0.3, 1.0/4.0)),\
                            1),\
        edges={'z10':["arc", arc_p1],\
        'z11':["arc", arc_p2],\
        'x00':["spline", (spline_p11, spline_p12)],\
        'x10':["spline", (spline_p21, spline_p22)]},\
        mergePatchPairs = [], metric='m'
        )
#print(b1)

bmd=BlockMeshDict(b1)
for i in range(72):
    # 1. geometric transformations, such as translation, rotation,
    #    reflection and scaling, are implemented for vertices, edges, faces,
    #    blocks and BlockMeshDicts
    # 2. Two or more neighboring BlockMeshDict objects can be merged by 
    #    adding them together
    bmd += rotate(b1,center=(0,0,0),axis=(0,0,1),angle= math.radians(wd)*i)
print(bmd)
