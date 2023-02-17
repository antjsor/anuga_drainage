import anuga, anuga.parallel, numpy, time, os, glob
from anuga.operators.rate_operators import Polygonal_rate_operator
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain, Inlet_operator
import anuga.utilities.spatialInputUtil as su
from anuga import distribute, myid, numprocs, finalize, barrier
from anuga.parallel.parallel_operator_factory import Inlet_operator, Boyd_box_operator, Boyd_pipe_operator
from anuga import Rate_operator
from anuga import Region


import numpy as np
from hymo import SWMMInpFile
import pickle
from pyswmm import Simulation, Nodes, Links
import matplotlib.pyplot as plt

import time
from coupling_functions.inlet_initialization import initialize_inlets
from coupling_functions.coupling import calculate_Q

time_average = 10 # sec
dt           = 1.0     # yield step
ft           = 25 # final timestep
# input_rate   = 0.102
input_rate   = 0.0


from pyswmm import SystemStats


output_frequency = 1
do_print         = True
do_data_save     = False

basename = 'model/terrain'
inp_name = 'real_example_inflow_1cms.inp'
outname  = 'real_example_swmm_1cms_in'
meshname = 'model/terrain.tsh'

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------
riverWall_csv_files = glob.glob('model/wall/*.csv') # Make a list of the csv files in BREAKLINES

(riverWalls, riverWall_parameters) = su.readListOfRiverWalls(riverWall_csv_files)

CatchmentDictionary = {'model/kerb/kerb1.csv':0.01, 'model/kerb/kerb2.csv':0.01}
    
bounding_polygon = anuga.read_polygon('model/domain.csv')
interior_regions = anuga.read_polygon_dir(CatchmentDictionary, 'model/kerb')


create_mesh_from_regions(bounding_polygon,
    boundary_tags={'inflow': [12], 'bottom': [0,1,2,3,4,5], 'top': [7,8,9,10,11], 'outflow': [6]},
    #boundary_tags=None,
    maximum_triangle_area = 0.1,
    breaklines = riverWalls.values(),
    interior_regions = interior_regions,
    filename = meshname,
    use_cache = False,
    verbose = False)

#------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------
domain = anuga.Domain(meshname, use_cache=False, verbose=False)
domain.set_minimum_storable_height(0.0)
domain.riverwallData.create_riverwalls(riverWalls,verbose = False) 
domain.set_name(outname) 
# 
#------------------------------------------------------------------------------
# APPLY MANNING'S ROUGHNESSES
#------------------------------------------------------------------------------

domain.set_quantity('friction', 0.025)
domain.set_quantity('stage', 0)
domain.set_quantity('elevation', filename=basename+'.csv', use_cache=False, verbose=False, alpha=0.99)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)  
Bd = anuga.Dirichlet_boundary([0,0,0])

domain.set_boundary({'inflow': Br, 'bottom': Br, 'outflow': Bd, 'top': Br})
 
# ------------------------------------------------------------------------------
# Setup inject water
# ------------------------------------------------------------------------------

input1_anuga_region   = Region(domain, radius=1.0, center=(305694.91,6188013.94))
input1_anuga_inlet_op = Inlet_operator(domain, input1_anuga_region, Q=input_rate) 

sim = Simulation(inp_name)
inp = SWMMInpFile(inp_name)


link_volume_0 = 0
for link in Links(sim):
    link_volume_0 += link.volume

# old_inlet_vol = [node.volume for node in Nodes(sim) if node.is_junction()]
node_ids      = [node.nodeid for node in Nodes(sim)]
in_node_ids   = [node.nodeid for node in Nodes(sim) if node.is_junction()]
n_in_nodes    = len(in_node_ids)


### Initialize inlet operators
inlet_area     = np.full((n_in_nodes),1.167)
Q_in_0         = n_in_nodes*[0.0]
Q_in_cumu      = sum(Q_in_0)
inlet_vol_cumu = 0
n_sides        = 6

inlet_operators,inlet_elevation,_,_ = initialize_inlets(domain,sim,inp,n_sides,inlet_area,Q_in_0,rotation = 0)

inlet_weir_length = 2*np.sqrt(np.pi*inlet_area)


Q_in_old       = np.zeros_like(inlet_elevation)
outfall_vol    = 0

times          = []
losses         = []

if do_data_save:
    Q_ins          = []
    conduit_depths = []
    node_heads     = []
    cumulative_inlet_flooding = np.array(n_in_nodes*[0.0])
    cumulative_inlet_flow     = np.array(n_in_nodes*[0.0])

system_routing = SystemStats(sim)

wall_clock_start = time.perf_counter()
sim.start()
old_inlet_vol = [- node.statistics['lateral_infow_vol'] + node.statistics['flooding_volume'] for node in Nodes(sim) if node.is_junction()]
node_volume    = sum(old_inlet_vol)

for t in domain.evolve(yieldstep=dt, finaltime=ft):
    anuga_depths = np.array([inlet_operators[in_id].inlet.get_average_depth() for in_id in in_node_ids])

    # Reset link volume at every iteration and sum volumes
    link_volume = 0
    for link in Links(sim):
        link_volume += link.volume

    if do_data_save:
        conduit_depths.append(np.array([link.depth for link in Links(sim)]))

    inlet_head_swmm   = np.array([node.head for node in Nodes(sim) if node.is_junction()])

    Q_in       = calculate_Q(inlet_head_swmm, anuga_depths, inlet_elevation, inlet_weir_length, inlet_area) # inputs between manual and auto checked to be the same 20/09
    Q_in       = ((time_average - dt)*Q_in_old + dt*Q_in)/time_average
    Q_in_old   = Q_in.copy()
    Q_in_cumu += sum(Q_in)
    if do_data_save:
        Q_ins.append(Q_in.copy())
        node_heads.append(inlet_head_swmm)

    # Simulate sewer with flow input
    for node, Qin in zip(Nodes(sim), Q_in): 
        node.generated_inflow(Qin)

    sim.step_advance(dt) 
    sim.next()

    ### Using flow methods methods
    inlet_flow_method_vol = np.array([-node.lateral_inflow + node.flooding for node in Nodes(sim) if node.is_junction()])*dt
    
    ### Compute inlet flow using volumes

    inlet_flood_vol = [node.statistics['flooding_volume'] for node in Nodes(sim) if node.is_junction()]
    inlet_lat_vol   = [node.statistics['lateral_infow_vol'] for node in Nodes(sim) if node.is_junction()]
    
    inlet_vol       = [- node.statistics['lateral_infow_vol'] + node.statistics['flooding_volume'] for node in Nodes(sim) if node.is_junction()]

    inlet_flow      = [(new_vol - old_vol)/dt for new_vol,old_vol in zip(inlet_vol,old_inlet_vol)]
    inlet_vol_cumu += sum(inlet_vol)
    old_inlet_vol   = inlet_vol

    # Compute statistics and append data
    inlet_idx = 0
    for node in Nodes(sim):
        if node.is_junction():
            inlet_operators[node.nodeid].set_Q(inlet_flow[inlet_idx])
            inlet_idx += 1
    
    outfall_vol += Links(sim)['Conduit_4'].flow*dt

    sewer_volume         = link_volume + node_volume
    domain_volume        = domain.get_water_volume()
    sewer_volume         = link_volume + node_volume
    boundary_flow        = domain.get_boundary_flux_integral()
    total_volume_correct = t*input_rate + boundary_flow + link_volume_0
    # domain_sewer_vol     = domain_volume + sewer_volume
    total_volume_real    = domain_volume + sewer_volume + outfall_vol
    # total_volume_outfall_vol    = domain_volume + sewer_volume + outfall_vol

    loss = total_volume_real - total_volume_correct
    losses.append(loss)

    if do_data_save:
        cumulative_inlet_flooding += np.array([node.flooding for node in Nodes(sim) if node.is_junction()])
        cumulative_inlet_flow     += np.array(inlet_flow)*dt
    times.append(t)

    if domain.yieldstep_counter%output_frequency == 0 and do_print:
        pass
        print('\nt = ',t)
        print(f'Q_in = {Q_in}')
        # print(f'Link4 flow: {Links(sim)["Conduit_4"].flow}')
        # print(f'outfall_vol     : {outfall_vol:.4f}')
        # print(f'Outflow_routing : {system_routing.routing_stats["outflow"]:.4f}')


        print(f'External_inflow      : {system_routing.routing_stats["external_inflow"]:.4f}')
        print(f'sum(inlet_lat_vol)   : {sum(inlet_lat_vol):.4f}')
        print(f'sum(inlet_flood_vol) : {sum(inlet_flood_vol):.4f}')
        print(f'-sum(inlet_vol)      : { sum(inlet_vol):.4f}')
        print(f'Q_in_cumu            : {Q_in_cumu:.4f}')
        print(f'external_flow - Q_in_cumu : {system_routing.routing_stats["external_inflow"] - Q_in_cumu:.4f}')
        # print(f'-sum(inlet_vol) + Q_in_cumu -: {-sum(inlet_vol) + Q_in_cumu- sum(inlet_flood_vol):.4f}')


sim.report()
sim.close()
wall_clock_end = time.perf_counter()
print('\n')
print(40*'#')
print(f'Loss : {loss}')
print(f'\nComputation time: {wall_clock_end - wall_clock_start:.2f} seconds')

# print(f'Loss = {loss:.2f}m^3 of total {t*input_rate}m^3')


if do_data_save:

    pick = outname + '.dat'

    data = {'times':times, 'conduit_depths':conduit_depths, 'node_heads':node_heads, 'Q_in':Q_in, 'inlet_vol_sum':sum(inlet_vol) }

    with open(pick, "wb") as f:
        pickle.dump(data, f)