import sys
sys.path.insert(0, '/home/llandsmeer/Projects/MEP/io/Glomerulus/arbormax/build/python')


from pathlib import Path

import arbor as arb
import matplotlib.pyplot as plt

class recipe(arb.recipe):
    def __init__(self):
        arb.recipe.__init__(self)
        self.tree = arb.segment_tree()
        # radius .5 and height 1 s.t. area is 4.71
        # on change, use wolfram alpha and fill in to get the right cm
        # 281 pF / 4.71 mm^2 = x F/m2
        self.tree.append(arb.mnpos, (0, 0, 0, .5), (1, 0, 0, .5), 1)
        self.props = arb.neuron_cable_properties()
        self.cat = arb.IOU_catalogue()
        self.props.register(self.cat)
        print(list(self.cat))
        d = arb.decor()
        d.paint('(all)', arb.mechanism('adex', dict(I=1.5e-9)))
        d.paint('(all)', cm=1)
        d.set_property(Vm=-60)
        labels = arb.label_dict()
        labels['root'] = '(root)'
        self.cell = arb.cable_cell(self.tree, labels, d)

    def global_properties(self, _):
        return self.props

    def num_cells(self):
        return 1

    def num_targets(self, gid):
        return 0

    def probes(self, gid):
        return [
                arb.cable_probe_membrane_voltage('"root"'),
                arb.cable_probe_density_state('"root"', "adex", "w"),
                ]

    def num_sources(self, gid):
        return 0

    def cell_kind(self, gid):
        return arb.cell_kind.cable

    def cell_description(self, gid):
        return self.cell

rcp = recipe()
ctx = arb.context()
dom = arb.partition_load_balance(rcp, ctx)
sim = arb.simulation(rcp, dom, ctx)

handle = sim.sample((0, 0), arb.regular_schedule(0.025))
handle_ = sim.sample((0, 1), arb.regular_schedule(0.025))
sim.run(dt=0.025, tfinal=200)
t, vs = sim.samples(handle)[0][0].T
t, w = sim.samples(handle_)[0][0].T
w *= 1e9
print(vs)
print(w)

plt.plot(t, vs, label='vm')
plt.plot(t, w, label='W')
plt.legend()
plt.show()
