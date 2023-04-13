


def calculate_Q(head1D, depth2D, bed2D, length_weir, area_manhole, cw=0.67, co=0.67, eps=1e-14):
    """
    Routine to calculate coupling discharge between 2D and 1D models

    Calculations based on weir and orifice equations.

    Reference:
    A methodology for linking 2D overland flow models with the sewer network model SWMM 5.1 
    based on dynamic link libraries
    Leandro, Jorge, Martins, Ricardo
    Water Science and Technology
    2016, 73, 3017-3026

    cw is the weir discharge coefficient, 
    length_weir is the weir crest width [m], 
    area_manhole is the manhole area [m2] 
    co is the orifice discharge coefficient.

    """

    import numpy as np
    from anuga import g

    with np.errstate(invalid='ignore'):
        Q = np.zeros_like(head1D)

        #print(Q)
        # if head1D < bed2D use Weir Equation (Reference Equation (10)):
        Q = np.where(head1D<bed2D, cw*length_weir*depth2D*np.sqrt(2*g*depth2D), Q)

        #print(Q)
        # If head1D > bed2D and  head1D < (depth2D + bed2D) use orifice equation (Equation (11))
        Q = np.where(np.logical_and(bed2D<=head1D, head1D<depth2D+bed2D) , co*area_manhole*np.sqrt(2*g*(depth2D+bed2D-head1D)), Q)

        #print(Q)

        #print(head1D,depth2D, bed2D)
        #print(head1D-depth2D-bed2D)
        
        # Otherwise if h1d >= (depth2D + bed2D) use orifice equation (Equation (11)) surcharge
        Q = np.where(head1D>depth2D+bed2D+eps,  -co*area_manhole*np.sqrt(2*g*(head1D-depth2D-bed2D)), Q)
        #print(Q)

    return Q

def calculate_Q_v2(head1D, depth2D, bed2D, length_weir, area_manhole, inlet_volumes, dt, cw=0.67, co=0.67, eps=1e-14):
    """
    Routine to calculate coupling discharge between 2D and 1D models

    Calculations based on weir and orifice equations.

    Reference:
    A methodology for linking 2D overland flow models with the sewer network model SWMM 5.1 
    based on dynamic link libraries
    Leandro, Jorge, Martins, Ricardo
    Water Science and Technology
    2016, 73, 3017-3026

    cw is the weir discharge coefficient, 
    length_weir is the weir crest width [m], 
    area_manhole is the manhole area [m2] 
    co is the orifice discharge coefficient.

    Modified to never take out more water than exists on terrain above the node.

    """

    import numpy as np
    from anuga import g

    with np.errstate(invalid='ignore'):
        Q = np.zeros_like(head1D)

        # # if head1D < bed2D use Weir Equation (Reference Equation (10)):
        # Q = np.where(head1D<bed2D, np.array([min(flow*dt,volume) for flow, volume in zip (cw*length_weir*depth2D*np.sqrt(2*g*depth2D),inlet_volumes)]), Q)
        
        # # If head1D > bed2D and  head1D < (depth2D + bed2D) use orifice equation (Equation (11))
        # Q = np.where(np.logical_and(bed2D<=head1D, head1D<depth2D+bed2D) , np.array([min(flow*dt,volume) for flow, volume in zip (co*area_manhole*np.sqrt(2*g*(depth2D+bed2D-head1D)),inlet_volumes)]), Q)

        # if head1D < bed2D use Weir Equation (Reference Equation (10)):
        Q = np.where(head1D<bed2D, np.array([min(flow,volume/dt) for flow, volume in zip (cw*length_weir*depth2D*np.sqrt(2*g*depth2D), inlet_volumes)]), Q)
        
        # If head1D > bed2D and  head1D < (depth2D + bed2D) use orifice equation (Equation (11))
        Q = np.where(np.logical_and(bed2D<=head1D, head1D<depth2D+bed2D) , np.array([min(flow,volume/dt) for flow, volume in zip (co*area_manhole*np.sqrt(2*g*(depth2D+bed2D-head1D)), inlet_volumes)]), Q)

        # Otherwise if h1d >= (depth2D + bed2D) use orifice equation (Equation (11)) surcharge
        Q = np.where(head1D>depth2D+bed2D+eps,  -co*area_manhole*np.sqrt(2*g*(head1D-depth2D-bed2D)), Q)
        #print(Q)

    return Q


# def node_link_connectivity(sim):
#     link2node = dict()
#     node2link = {node.nodeid: [] for node in Nodes(sim)}
#     for link in Links(sim):
#         link2node[link.linkid] = link.connections
#         for nodeid in link.connections:
#             node2link[nodeid].append(link.linkid)
#     return link2node, node2link

# def inlet_lateral_inflow(sim):
#     link2node, node2link = node_link_connectivity(sim)
    