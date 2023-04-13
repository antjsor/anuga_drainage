import math
from pyswmm import Nodes
from hymo import SWMMInpFile
from anuga import Inlet_operator,Region
import numpy as np

def n_sided_inlet(n_sides, area, inlet_coordinate, rotation):
    # Computes the vertex coordinates and side length of a regular polygon with:
    # Number of sides = n_sides
    # Area = area
    #TODO Would it be more sensible to smaple points along circle and just use paramters for desired circle areas in weir/orifice equations?
    if n_sides < 3:
        raise RuntimeError('A polygon should have at least 3 sides')

    one_segment = math.pi * 2 / n_sides
    side_length = math.sqrt(4.0*area*math.tan(math.pi/n_sides)/n_sides)
    
    radius = side_length/(2.0*math.sin(math.pi/n_sides))

    vertex = [
        (math.sin(one_segment * i + rotation) * radius,
        math.cos(one_segment * i + rotation) * radius)
        for i in range(n_sides)]

    vertex = [[sum(pair) for pair in zip(point, inlet_coordinate)]
            for point in vertex]
        
    return vertex, side_length




def initialize_inlets(domain, inp, inlet_ids = None, n_sides = 6, manhole_areas = [1], Q_in_0 = 0.0, rotation = -np.pi/4, zero_velocity = False, expand_polygon = True):
    if n_sides < 3:
        raise RuntimeError('A polygon should have at least 3 sides')

    if not(isinstance(manhole_areas, int) or isinstance(manhole_areas,float) or isinstance(manhole_areas,list) or isinstance(manhole_areas,np.ndarray)):
        raise RuntimeError('Invalid masnhole area format')
    if not(inlet_ids):
        raise RuntimeError('Must supply at least one node id for inlet operator')

    inlet_operators = dict()
    elevation_list  = []
    circumferences  = []
    polygons        = []

    ## TODO replace by and statement
    if isinstance(manhole_areas, list) or isinstance(manhole_areas, np.ndarray):
        if len(manhole_areas) != len(inlet_ids):
            raise RuntimeError('Number of areas must be equal to number of node ids or a single float value ')
        else:
            pass
    if isinstance(Q_in_0, list) or isinstance(Q_in_0, np.ndarray):
        if len(Q_in_0) != len(inlet_ids):
            raise RuntimeError('Number of inflows must be equal to number of node ids or single float value ')
        else:
            pass

    elif not(isinstance(manhole_areas, int) or isinstance(manhole_areas,float)):
        raise RuntimeError('Invalid manhole area data type')

    for inlet_idx, nodeid in enumerate(inlet_ids):
        inlet_coordinates = [inp.coordinates.loc[nodeid].X_Coord, inp.coordinates.loc[nodeid].Y_Coord]
        
        if isinstance(manhole_areas ,list) or isinstance(manhole_areas, np.ndarray):
            vertices, side_length = n_sided_inlet(n_sides, manhole_areas[inlet_idx], inlet_coordinates, rotation)
        else:
            vertices, side_length = n_sided_inlet(n_sides, manhole_areas, inlet_coordinates, rotation)
    
        if isinstance(Q_in_0 ,list) or isinstance(Q_in_0, np.ndarray):
            inlet_operators[nodeid] = Inlet_operator(domain, Region(domain,polygon = vertices,expand_polygon = expand_polygon), Q_in_0[inlet_idx], zero_velocity=zero_velocity)
        else:
            inlet_operators[nodeid] = Inlet_operator(domain, Region(domain,polygon = vertices,expand_polygon = expand_polygon), Q_in_0, zero_velocity=zero_velocity)


        elevation_list.append(inlet_operators[nodeid].inlet.get_average_elevation())
        polygons.append(vertices)
        circumferences.append(side_length)

    polygons         = np.array(polygons)
    inlet_elevations = np.array(elevation_list)
    circumferences   = np.array(circumferences)

    return inlet_operators, inlet_elevations, circumferences, polygons

# def initialize_inlets(domain, sim, inp, n_sides = 6, manhole_areas = [1],
#                        Q_in_0 = [1], rotation = -np.pi/4, zero_velocity = True, expand_polygon = True):
#     print(f'zero_velocity: {zero_velocity}')
#     print(f'expand_polygon: {expand_polygon}')
#     print(domain)

#     if n_sides < 3:
#         raise RuntimeError('A polygon should have at least 3 sides')

#     if not(isinstance(manhole_areas, int) or isinstance(manhole_areas,float) or isinstance(manhole_areas,list) or isinstance(manhole_areas,np.ndarray)):
#         raise RuntimeError('Invalid ')

#     inlet_operators = dict()
#     elevation_list  = []
#     circumferences  = []
#     polygons        = []
#     in_nodes        = [node for node in Nodes(sim) if node.is_junction()]

#     for inlet_idx, node in enumerate(in_nodes):
#         if isinstance(manhole_areas,list) or isinstance(manhole_areas,np.ndarray):
#             inlet_area = manhole_areas[inlet_idx] 

#         inlet_coordinates = [inp.coordinates.loc[node.nodeid].X_Coord, inp.coordinates.loc[node.nodeid].Y_Coord]
#         vertices, side_length = n_sided_inlet(n_sides, inlet_area, inlet_coordinates, rotation)
        
#         inlet_operators[node.nodeid] = Inlet_operator(domain, Region(domain,polygon = vertices, expand_polygon = expand_polygon), Q_in_0[inlet_idx], zero_velocity=zero_velocity)

#         elevation_list.append(inlet_operators[node.nodeid].inlet.get_average_elevation())
#         circumferences.append(n_sides*side_length)
#         polygons.append(vertices)

#     polygons         = np.array(polygons)
#     inlet_elevations = np.array(elevation_list)
#     circumferences   = np.array(circumferences)

#     return inlet_operators, inlet_elevations, circumferences, polygons

