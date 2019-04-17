# for compatibility to Py2.7
from __future__ import unicode_literals, print_function
from six import string_types
import io, re
from collections import Iterable, Counter
from copy import copy, deepcopy
from string import Template
from .euclid import Vector3
from itertools import groupby

# Vertex is subclass of Vector3 from the pyeuclid library 
class Vertex(Vector3):
    def __init__(self, x, y=None, z=None):
        if (y is not None) and (z is not None):
            self.x = x
            self.y = y
            self.z = z
        else:
            self.x = x[0]
            self.y = x[1]
            self.z = x[2]

    def __repr__(self):
        return '(%18.15g %18.15g %18.15g)' % self._grid()

    def _grid(self, coords = None):
        # Grid resolution (GR) in current metric
        GR = 1e-9
        # Any geometry details that are smaller than GR will be lost.
        # With a reasonable choice of metric, GR = 1e-9 should be
        # enough for general problems.
        if coords is None:
            return (round(self.x / GR) * GR,\
                    round(self.y / GR) * GR,\
                    round(self.z / GR) * GR)
        else:
            return (round(coords[0] / GR) * GR,\
                    round(coords[1] / GR) * GR,\
                    round(coords[2] / GR) * GR)

    def __hash__(self):
        return hash(self._grid())

    def __eq__(self, other):
        assert hasattr(other, '__len__') and len(other) == 3
        return self._grid() == \
                self._grid((other[0], other[1], other[2]))
    def _transform(self, method, *args, **kwargs):
        if method == "rotate":
            center = Vertex(kwargs['center'])
            axis = Vertex(kwargs['axis'])
            angle = kwargs['angle']
            geometry_moveToOrigin = self - center
            geometry_applyRotate = \
                    geometry_moveToOrigin.rotate_around(axis, angle)
            geometry_final = geometry_applyRotate + center
        if method == "scale":
            center = Vertex(kwargs['center'])
            ratio = kwargs['ratio']
            geometry_final = (self - center) * ratio + center
        if method == "mirror":
            center = Vertex(kwargs['center'])
            normal = Vertex(kwargs['normal'])
            geometry_final = (self - center).reflect(normal) + center
        if method == "translate":
            vector = Vertex(kwargs['vector'])
            geometry_final = self + vector
        return  Vertex(\
                geometry_final[0], \
                geometry_final[1], \
                geometry_final[2]
                )

class Edge(object):
    def __init__(self, start=(0,0,0), end=(0,0,0), edgeType="line", \
            points=[], index=None):
        if len(points)>0 and \
                (isinstance(points[0],float) or isinstance(points[0],int)):
            self.vertices = [tuple(start), tuple(end), tuple(points)]
        else:
            self.vertices = [tuple(start), tuple(end)]
            for p in points:
                self.vertices.append(tuple(p))
        self._type = edgeType

    def __repr__(self):
        start = Vertex(self.vertices[0])
        end = Vertex(self.vertices[1])
        pList = self.vertices[2:]
        if len(pList) == 1:
            points = Vertex(pList[0])
        else:
            points = foamList([Vertex(p) for p in pList])
        return '%s %s %s %s' %\
                (self._type,str(start),str(end),str(points))

    def __copy__(self):
        e = Edge()
        e.vertices = [copy(v) for v in self.vertices]
        e._type = self._type
        return e
    copy = __copy__

    def __hash__(self):
        start = Vertex(self.vertices[0])._grid()
        end = Vertex(self.vertices[1])._grid()
        return hash(tuple(sorted([start,end])))

    def __eq__(self, other):
        assert isinstance(other, Edge)
        s1 = Vertex(self.vertices[0])._grid()
        s2 = Vertex(self.vertices[1])._grid()
        o1 = Vertex(other.vertices[0])._grid()
        o2 = Vertex(other.vertices[1])._grid()
        return tuple(sorted([s1, s2])) == tuple(sorted([o1, o2]))

    def _transform(self, method, *args, **kwargs):
        newEdge = copy(self)
        newEdge.vertices = [transform(Vertex(v), method, *args, **kwargs) \
                for v in self.vertices]
        return newEdge

class Face(object):
    def __init__(self, vertices=[Vertex(0,0,0),]*4, \
            name="empty_defaultFaces", index=None):
        self.vertices = vertices
        # Personally, I like the "patchType_patchName" naming convention,
        # for example, a wall patch XXXX could be named as "wall_XXXX"
        self._name = '_'.join(name.split('_')[1:])
        self._type = name.split('_')[0]
        assert self._type in ['patch', 'wall', 'symmetryPlane', 'symmetry',\
                'empty', 'wedge', 'cyclic', 'cyclicAMI', 'processor']

    def __repr__(self):
        out = "( "
        for v in self.vertices:
            out += str(Vertex(v)) + " "
        out += ")"
        return out

    def __copy__(self):
        f = Face()
        f.vertices = [copy(v) for v in self.vertices]
        f._name = self._name
        f._type = self._type
        return f
    copy = __copy__

    def __hash__(self):
        return  hash(tuple(sorted(\
                [Vertex(v)._grid() for v in self.vertices]
                )))

    def __eq__(self, other):
        assert  isinstance(other, Face)
        return  tuple(sorted(\
                [Vertex(v)._grid() for v in self.vertices]))\
                ==\
                tuple(sorted(\
                [Vertex(v)._grid() for v in other.vertices]))

    def _transform(self, method, *args, **kwargs):
        newFace = copy(self)
        newFace.vertices = [transform(Vertex(v), method, *args, **kwargs) \
                for v in self.vertices]
        return newFace

class Patch(object):
    def __init__(self, faces):
        self.faces = faces
        nameSet = set([f._name for f in faces])
        typeSet = set([f._type for f in faces])
        # All faces must have the same name and type
        assert len(nameSet) == 1
        self._name = list(nameSet)[0]
        assert len(typeSet) == 1
        self._type = list(typeSet)[0]

    def __repr__(self):
        out = self._name + "\n    {\n        type " + self._type + \
                ";\n        faces\n        (\n" 
        for f in self.faces:
            out += "            " + str(f) + "\n"
        out += "        );\n    }"
        return out

class MergePatchPair(object):
    def __init__(self, mergePatchPair=None):
        self._mergepatchpair = \
                tuple(str(mergePatchPair[0]),str(mergePatchPair[1]))

    def __repr__(self):
        if self._mergepatchpair is not None:
            return "(" + self._mergepatchpair[0] + " "\
                    + self._mergepatchpair[1] + ")"
        else:
            return ""

    def __copy__(self):
        mpp = MergePatchPair()
        mpp._mergepatchpair = copy(self._mergepatchpair)
        return mpp
    copy = __copy__

class Block(object):
    def __init__(self, var, region=None, cells=None, \
            grading = None, edges=None, boundaries=None):
        # Prepare args
        passed_args = {'region':region, 'cells':cells, \
            'grading':grading, 'edges':edges, 'boundaries':boundaries}
        default_args = {'region':"fluid", 'cells':(1,1,1), \
            'grading':SimpleGrading(1,1,1), 'edges':{}, 'boundaries':{}}
        
        # Build vertices
        if isinstance(var, Block):
            self.vertices = [Vertex(v) for v in var.vertices]
            copied_args = {'region':var.region, 'cells':var.cells, \
                'grading':var.grading, 'edges':var.edges, \
                'boundaries':var.faces}
        else:
            assert len(var) == 8
            for v in var:
                assert isinstance(v, Vertex) or len(v) == 3
            self.vertices = [Vertex(v) for v in var]
            copied_args = {'region':None, 'cells':None, 'grading':None, \
                'edges':None, 'boundaries':None}

        # Initialize variables, use passed args first, then try copied args
        # finally fallback to the default args

        # Region
        if passed_args['region'] is not None:
            self.region = passed_args['region']
        else:
            if copied_args['region'] is not None:
                self.region = copied_args['region']
            else:
                self.region = default_args['region']
        
        # Cells 
        if passed_args['cells'] is not None:
            self.cells = passed_args['cells']
        else:
            if copied_args['cells'] is not None:
                self.cells = copied_args['cells']
            else:
                self.cells = default_args['cells']

        # Grading
        if passed_args['grading'] is not None:
            self.grading = deepcopy(passed_args['grading'])
        else:
            if copied_args['grading'] is not None:
                self.grading = deepcopy(copied_args['grading'])
            else:
                self.grading = deepcopy(default_args['grading'])

        # Edges
        # 1. initialization 
        edge_args = {\
            'x00':{'index':(0, 1)}, 'x10':{'index':(3, 2)},\
            'x11':{'index':(7, 6)}, 'x01':{'index':(4, 5)},\
            'y00':{'index':(0, 3)}, 'y10':{'index':(1, 2)},\
            'y11':{'index':(5, 6)}, 'y01':{'index':(4, 7)},\
            'z00':{'index':(0, 4)}, 'z10':{'index':(1, 5)},\
            'z11':{'index':(2, 6)}, 'z01':{'index':(3, 7)}}
        for k, v in edge_args.items():
            v['start'] = self.vertices[v['index'][0]]
            v['end']   = self.vertices[v['index'][1]]
        # 2. generate edges
        if passed_args['edges'] is not None:
            edges = passed_args['edges']
            for k, v in edges.items():
                edge_args[k]['edgeType'] = v[0]
                edge_args[k]['points'] = v[1]
            self.edges = {k : Edge(**args) for k, args in edge_args.items()}
        else:
            if copied_args['edges'] is not None:
                self.edges = {k : copy(v) for k, v in var.edges.items()}
            else:
                self.edges = {k : Edge(**args) \
                    for k, args in edge_args.items()}

        # Boundaries
        # 1. initialization
        face_args = {\
            'xm':{'index':(0, 4, 7, 3)}, 'xp':{'index':(1, 2, 6, 5)},\
            'ym':{'index':(0, 1, 5, 4)}, 'yp':{'index':(2, 3, 7, 6)},\
            'zm':{'index':(0, 3, 2, 1)}, 'zp':{'index':(4, 5, 6, 7)}} 
        for k, v in face_args.items():
            v['vertices'] = tuple([self.vertices[v['index'][i]] \
                    for i in range(4)])
        # 2. generate faces
        if passed_args['boundaries'] is not None:
            faces = passed_args['boundaries']
            for k, v in faces.items():
                face_args[k]['name'] = v
                self.faces = {k : Face(**args) \
                    for k, args in face_args.items()}
        else:
            if copied_args['boundaries'] is not None:
                self.faces={k : copy(v) for k, v in var.faces.items()}
            else:
                self.faces = {k : Face(**args) \
                    for k, args in face_args.items()}

    def __repr__(self):
        out = "hex ( "
        for v in self.vertices:
            out += str(Vertex(v))+" "
        out += ")"
        out += " " + self.region
        out += " (" + str(self.cells[0]) + " " + str(self.cells[1])\
                + " " + str(self.cells[2]) + ") " 
        out += self.grading.format()
        return out

    def __hash__(self):
        return  hash(tuple(sorted(\
                [Vertex(v)._grid() for v in self.vertices]
                )))

    def __eq__(self, other):
        assert  isinstance(other, Block)
        return  tuple(sorted(\
                [Vertex(v)._grid() for v in self.vertices]))\
                ==\
                tuple(sorted(\
                [Vertex(v)._grid() for v in other.vertices]))

    def __copy__(self):
        b = Block([copy(v) for v in self.vertices])
        b.region = self.region
        b.cells = tuple(self.cells)
        b.grading = deepcopy(self.grading)
        b.faces = {k : copy(face) for k, face in self.faces.items()}
        b.edges = {k : copy(edge) for k, edge in self.edges.items()}
        return b
    copy = __copy__

    def _transform(self, method, *args, **kwargs):
        newBlock = copy(self)
        newBlock.vertices = [transform(Vertex(v), method, *args, **kwargs)\
                for v in self.vertices]
        newBlock.edges = {k : transform(e, method, *args, **kwargs) \
                for k, e in self.edges.items()}
        newBlock.faces = {k : transform(f, method, *args, **kwargs) \
                for k, f in self.faces.items()}
        if method == "mirror":
            # Adjust the vertices order to ensure the right-hand rule
            nv = tuple(newBlock.vertices)
            newBlock.vertices = [copy(nv[0]), copy(nv[3]), \
                    copy(nv[2]), copy(nv[1]),\
                    copy(nv[4]), copy(nv[7]),\
                    copy(nv[6]), copy(nv[5])]
        return newBlock

class BlockMeshDict(list):
    def __init__(self, var, region=None, cells=None, \
            grading = None, edges=None, boundaries=None, \
            metric = None, mergePatchPairs = None): 
        super(BlockMeshDict, self).__init__()

        # 1. Case 1, var is BlockMeshDict
        if isinstance(var, BlockMeshDict):
            # 1.1 set metric
            if metric is not None:
                self.set_metric(metric)
            else:
                self.convert_to_meters = var.convert_to_meters
            # 1.2 set mergePatchPairs
            if mergePatchPairs is not None:
                if isinstance(mergePatchPairs, list):
                    self.mergePatchPairs = mergePatchPairs
                else:
                    self.mergePatchPairs = [mergePatchPairs,]
            else:
                self.mergePatchPairs = var.mergePatchPairs
            # 1.3 set blocks 
            for b in var:
                self.append(b)
        else:
            # 2. Case 2, var is a block
            # 3. Case 3, var is a list of blocks
            # 4. Case 4, var is a list of vertices
            # *.1 set metric
            if metric is not None:
                self.set_metric(metric)
            else:
                self.set_metric('m')
            # *.2 set mergePatchPairs
            if mergePatchPairs is not None:
                if isinstance(mergePatchPairs, list):
                    self.mergePatchPairs = mergePatchPairs
                else:
                    self.mergePatchPairs = [mergePatchPairs,]
            else:
                self.mergePatchPairs = []
            # 2.3 Case 2, var is a block
            if isinstance(var, Block):
                self.append(Block(var, region=region, cells=cells, \
                        grading = grading, edges=edges, \
                        boundaries=boundaries))
            else:
                # 3.3 Case 3, var is a list of blocks
                if isinstance(var, list):
                    vertices = []
                    for v in var:
                        if isinstance(v, Block):
                            self.append(Block(v,region=region,cells=cells,\
                                    grading = grading, edges=edges, \
                                    boundaries=boundaries))
                        else:
                            # 4.3 Case 4, var is a list of vertices
                            if isinstance(v, Vertex):
                                vertices.append(v)
                    if len(set(vertices)) > 3: 
                    # it requires at least 4 vertices to build a 3D geometry
                        self.append(Block(vertices, region=region, \
                                cells=cells, grading=grading, \
                                edges=edges, boundaries=boundaries))

    def set_metric(self, metric):
        """set self.comvert_to_meters by word"""
        metricsym_to_conversion = {
            'km': 1000,
            'm': 1,
            'cm': 0.01,
            'mm': 0.001,
            'um': 1e-6,
            'nm': 1e-9,
            'A': 1e-10,
            'Angstrom': 1e-10}
        self.convert_to_meters = metricsym_to_conversion[metric]

    def add_mergePatchPairs(self, m):
        self.mergePatchPairs.append(mergePatchPair(m))

    def remove_mergePatchPairs(self, m):
        self.mergePatchPairs.remove(mergePatchPair(m))

    def build_geometry(self):
        vertices = []
        edges = []
        faces = []
        patches = []
        blocks = list(set(self)) # Remove duplicated blocks
        mergePatchPairs = list(set(self.mergePatchPairs))
        for block in blocks:
            vertices += block.vertices
            edges += filter(lambda e : e._type != "line", \
                    block.edges.values())
            faces += block.faces.values()
        # 1. Remove internal faces
        boundaryfaces = [face \
                for face, count in Counter(faces).items() if count<2]
        # 2. Build patch dictionary from face list
        sortedfaces = sorted(list(set(boundaryfaces)),\
                key = lambda f : f._name)
        for patch,faceNames in groupby(sortedfaces, lambda f:hash(f._name)):
            patches.append(Patch(list(faceNames)))
        return  foamList(set(vertices)), \
                foamList(set(edges)), \
                foamList(set(patches)), \
                foamList(blocks), \
                foamList(mergePatchPairs)
    
    def _transform(self, method, *args, **kwargs):
        newBMD = copy(self)
        newBMD.clear()
        for b in self:
            newBMD.append(transform(b, method, *args, **kwargs))
        return newBMD

    def __copy__(self):
        bmd = BlockMeshDict([Vertex(0,0,0),]*8)
        bmd.convert_to_meters = self.convert_to_meters
        bmd.mergePatchPairs = [copy(m) for m in self.mergePatchPairs]
        for b in self:
            bmd.append(b.copy())
        return bmd
    copy = __copy__

    def __add__(self, other):
        assert isinstance(other, BlockMeshDict)

        # We choose to warn the user if different metrics were used,
        assert self.convert_to_meters == other.convert_to_meters
        # alternatively, we can scale the added BlockMeshDict by
        # other = scale(other, center=(0,0,0), \
        #       ratio=other.convert_to_meters/self.convert_to_meters

        bmd = copy(self)
        for b in other:
            bmd.append(b)
        for m in other.mergePatchPairs:
            bmd.mergePatchPairs.append(m)
        return bmd

    def _iadd_(self, other):
        assert isinstance(other, BlockMeshDict)
        # See the comments above
        assert self.convert_to_meters == other.convert_to_meters
        for b in other:
            self.append(b)
        for m in other.mergePatchPairs:
            self.mergePatchPairs.append(m)

    def __repr__(self):
        vertices, edges, patches, blocks, mergePatchPairs =\
                self.build_geometry()
        # Build replacement dict: from vertex coords (x, y, z) to index 
        replacements = {str(k):str(v) for v, k in enumerate(vertices)}
        template = Template(r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  6.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters $metricconvert;

vertices
$vertices;

edges
$edges;

blocks
$blocks;

boundary
$patches;

mergePatchPairs
$mergepatchpairs;

// ************************************************************************* //
''')
        return  template.substitute(\
                metricconvert=str(self.convert_to_meters),\
                vertices=str(vertices),\
                edges=multireplace(str(edges), replacements),\
                blocks=multireplace(str(blocks), replacements),\
                patches=multireplace(str(patches), replacements),\
                mergepatchpairs=\
                multireplace(str(mergePatchPairs), replacements))

#  Customize the built-in list type
class foamList(list):
    def __init__(self, *args, **kwargs):
        super(foamList, self).__init__(*args, **kwargs)
    
    def __repr__(self):
        out = "(\n"
        for v in self:
            out += "    "+str(v)+"\n"
        out +=")"
        return out

# Geometry transformation
def transform(geometry, method, *args, **kwargs):
    return geometry._transform(method, *args, **kwargs)

# rotation
# args: center=(0,0,0), axis=(0,0,1), angle=0.0
def rotate(geometry, *args, **kwargs):
    return transform(geometry, "rotate", *args, **kwargs)

# scale
# args: center=(0,0,0), ratio=1
def scale(geometry, *args, **kwargs):
    return transform(geometry, "scale", *args, **kwargs)

# mirror
# args: center=(0,0,0), normal=(0,0,1)
def mirror(geometry, *args, **kwargs):
    return transform(geometry, "mirror", *args, **kwargs)

# translate
# args: vector=(0,0,0)
def translate(geometry, *args, **kwargs):
    return transform(geometry, "translate", *args, **kwargs)

# Python string multi-replacement code by bgusach from:
# https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
def multireplace(string, replacements):
    """
    Given a string and a replacement map, it returns the replaced string.
    :param str string: string to execute replacements on
    :param dict replacements: replacement dictionary {value to find: value to replace}
    :rtype: str
    """
    # Place longer ones first to keep shorter substrings from matching where the longer ones should take place
    # For instance given the replacements {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should produce
    # 'hey ABC' and not 'hey ABc'
    substrs = sorted(replacements, key=len, reverse=True)

    # Create a big OR regex that matches any of the substrings to replace
    regexp = re.compile('|'.join(map(re.escape, substrs)))

    # For each match, look up the new string in the replacements
    return regexp.sub(lambda match: replacements[match.group(0)], string)

# Grading classes from ofblockmeshdicthelper by takaakiaoki
# https://github.com/takaakiaoki/ofblockmeshdicthelper
class Grading(object):
    """base class for Simple- and Edge- Grading"""
    pass


class SimpleGradingElement(object):
    """x, y or z Element of simpleGrading. adopted to multi-grading
    """
    def __init__(self, d):
        """initialization
        d is single number for expansion ratio
          or iterative object consits (dirction ratio, cell ratio, expansion ratio)
        """
        self.d = d

    def format(self):
        if isinstance(self.d, Iterable):
            s = io.StringIO()
            s.write('( ')
            for e in self.d:
                s.write('( {0:g} {1:g} {2:g} ) '.format(e[0], e[1], e[2]))
            s.write(')')
            return s.getvalue()
        else:
            return str(self.d)


class SimpleGrading(Grading):
    """configutation for 'simpleGrading'
    """
    def __init__(self, x, y, z):
        if not isinstance(x, SimpleGradingElement):
            self.x = SimpleGradingElement(x)
        else:
            self.x = x
        if not isinstance(y, SimpleGradingElement):
            self.y = SimpleGradingElement(y)
        else:
            self.y = y
        if not isinstance(z, SimpleGradingElement):
            self.z = SimpleGradingElement(z)
        else:
            self.z = z

    def format(self):
        return 'simpleGrading ({0:s} {1:s} {2:s})'.format(self.x.format(), self.y.format(), self.z.format())

class EdgeGrading(Grading):
    """configutation for 'edgeGrading'
    """
    def __init__(self, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4):
        if not isinstance(x1, SimpleGradingElement):
            self.x1 = SimpleGradingElement(x1)
        else:
            self.x1 = x1
        if not isinstance(x2, SimpleGradingElement):
            self.x2 = SimpleGradingElement(x2)
        else:
            self.x2 = x2
        if not isinstance(x3, SimpleGradingElement):
            self.x3 = SimpleGradingElement(x3)
        else:
            self.x3 = x3
        if not isinstance(x4, SimpleGradingElement):
            self.x4 = SimpleGradingElement(x4)
        else:
            self.x4 = x4
        if not isinstance(y1, SimpleGradingElement):
            self.y1 = SimpleGradingElement(y1)
        else:
            self.y1 = y1
        if not isinstance(y2, SimpleGradingElement):
            self.y2 = SimpleGradingElement(y2)
        else:
            self.y2 = y2
        if not isinstance(y3, SimpleGradingElement):
            self.y3 = SimpleGradingElement(y3)
        else:
            self.y3 = y3
        if not isinstance(y4, SimpleGradingElement):
            self.y4 = SimpleGradingElement(y4)
        else:
            self.y4 = y4
        if not isinstance(x1, SimpleGradingElement):
            self.z1 = SimpleGradingElement(z1)
        else:
            self.z1 = z1
        if not isinstance(z2, SimpleGradingElement):
            self.z2 = SimpleGradingElement(z2)
        else:
            self.z2 = z2
        if not isinstance(z3, SimpleGradingElement):
            self.z3 = SimpleGradingElement(z3)
        else:
            self.z3 = z3
        if not isinstance(z4, SimpleGradingElement):
            self.z4 = SimpleGradingElement(z4)
        else:
            self.z4 = z4
        

    def format(self):
        return 'edgeGrading '\
                '({0:s} {1:s} {2:s} {3:s} '\
                '{4:s} {5:s} {6:s} {7:s} '\
                '{8:s} {9:s} {10:s} {11:s})'.format(
                    self.x1.format(), self.x2.format(), self.x3.format(), self.x4.format(),
                    self.y1.format(), self.y2.format(), self.y3.format(), self.y4.format(),
                    self.z1.format(), self.z2.format(), self.z3.format(), self.z4.format()
                    )
