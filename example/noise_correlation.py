#!/usr/bin/env python3

import arbor
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

def make_cable_cell(gid):
    tree = arbor.segment_tree()
    s = tree.append(arbor.mnpos, arbor.mpoint(-12, 0, 0, 6), arbor.mpoint(0, 0, 0, 6), tag=1)
    b0 = tree.append(s, arbor.mpoint(0, 0, 0, 2), arbor.mpoint(50, 0, 0, 2), tag=3)
    b1 = tree.append(b0, arbor.mpoint(50, 0, 0, 2), arbor.mpoint(50+50/sqrt(2), 50/sqrt(2), 0, 0.5), tag=3)
    b2 = tree.append(b0, arbor.mpoint(50, 0, 0, 1), arbor.mpoint(50+50/sqrt(2), -50/sqrt(2), 0, 1), tag=3)

    labels = arbor.label_dict()
    labels['soma'] = '(tag 1)'
    labels['dend'] = '(tag 3)'

    labels['synapse_site'] = '(location 1 0.5)'
    labels['root'] = '(root)'

    decor = arbor.decor()

    decor.paint('"soma"', 'pas')
    decor.paint('"dend"', 'pas')

    decor.place('"gjpos"', arbor.gap_junction(), 'a')
    labels['gjpos'] = '(root)'

    decor.paint('"soma"', arbor.mechanism('ou_noise/theta=1,alpha=0.5,sigma=1'))

    #decor.place('"synapse_site"', 'expsyn', 'syn')
    #decor.place('"root"', arbor.spike_detector(-10), 'detector')

    cell = arbor.cable_cell(tree, labels, decor)

    return cell

class ring_recipe (arbor.recipe):

    def __init__(self, ncells):
        arbor.recipe.__init__(self)
        self.ncells = ncells
        self.props = arbor.neuron_cable_properties()
        self.cat = arbor.default_catalogue()
        self.cat.extend(arbor.IOU_catalogue(), '')
        self.props.register(self.cat)

    def num_cells(self): return self.ncells
    def cell_description(self, gid): return make_cable_cell(gid)
    def cell_kind(self, gid): return arbor.cell_kind.cable
    def connections_on(self, gid): return []
    def num_gap_junction_sites(self, gid): return self.ncells-1
    def event_generators(self, gid): return []
    def probes(self, gid): return [arbor.cable_probe_membrane_voltage('"root"')]
    def global_properties(self, kind): return self.props
    def gap_junctions_on(self, gid):
        conns = []
        for i in range(self.ncells):
            if i == gid: continue
            conn = arbor.gap_junction_connection(local='a', peer=(i, 'a'), ggap=0.0)
            conns.append(conn)
        return conns



ncells = 128
recipe = ring_recipe(ncells)
context = arbor.context()
#context = arbor.context(gpu_id=0)
decomp = arbor.partition_load_balance(recipe, context)
sim = arbor.simulation(recipe, decomp, context)
sim.record(arbor.spike_recording.all)
handles = [sim.sample((gid, 0), arbor.regular_schedule(0.1)) for gid in range(ncells)]
sim.run(100, dt=0.01)

vs = []
for gid in range(ncells):
    samples, meta = sim.samples(handles[gid])[0]
    #plt.plot(samples[:,0], samples[:,1])
    m = np.isnan(samples[:,1])
    vs.append(samples[len(samples)//2:,1])
vs = np.array(vs)
vs = (vs - vs.mean(1)[:,None]) / vs.std(1)[:,None]
print(vs.shape)
img = abs(np.cov(vs))

plt.imshow(img)
plt.colorbar()

plt.show()

for v in vs:
    plt.plot(v, color='black', alpha=0.1)

plt.axis('off')
plt.tight_layout()

plt.show()
