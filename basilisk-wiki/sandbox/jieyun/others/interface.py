import numpy as np
"""
# Transform facets output into Tecplot format
"""

class Facets():
    def __init__(self) -> None:
        self.xy = [] # vertex coordinates
        self.nv = [] # number of vertices of each facet

    def get_data_facets(self, path):
        """ Read facets from file output by the output_facets function in Basilisk. """
        xy = []
        self.nv = []
        with open(path, mode='r') as f:
            nv = 0
            for line in f:
                if line:
                    xx = [float(tmp) for tmp in line.split()]
                    if (len(xx) > 0):
                        nv += 1
                        xy.append(xx)
                    else:
                        # finish reading one facet
                        self.nv.append(nv)
                        nv = 0

        self.xy = np.array(xy)
        print("Read facets from file:{:s}, number of points:{:d}, \
              number of elements:{:d}".format(path, len(self.xy), len(self.nv)))

    def show_facets(self, ax):
        nv = 0
        xy = self.xy
        for i in self.nv:
            _xy = np.zeros((i + 1, 3))
            _xy[:i,:] = xy[nv:nv+i, :]
            _xy[i, :] = _xy[0, :]
            ax.plot(_xy[:, 0], _xy[:, 1], _xy[:, 2],'r-', linewidth=1)
            nv += i
        ax.set_box_aspect((np.ptp(xy[:, 0]), np.ptp(xy[:, 1]), np.ptp(xy[:, 2])))

    def output_tecplot(self, path):
        """ Triangulate facets and output into tecplot format for visualization. """
        xy = self.xy
        _xy = []
        _edge = []  # connectivity information
        nv = 0
        nvertex = 0
        nelement = 0

        for i in self.nv:
            for xx in xy[nv:nv+i]:
                _xy.append(xx)
            xc = np.mean(xy[nv:nv+i,:], axis=0)
            if i > 3:
                # decompose a polygon with N edges into N triangles.
                _xy.append(xc)
                for ii in range(i):
                    _edge.append([nvertex + ii + 1, nvertex + (ii + 1) % i + 1, nvertex + i + 1])
                nvertex += (i + 1)
                nelement += i
            else:
                _edge.append([nvertex + 1, nvertex + 2, nvertex + 3])
                nvertex += i
                nelement += 1

            nv += i

        with open(path, mode='w') as f:
            f.write('TITLE = "NONE"\n')
            f.write('VARIABLES = "X", "Y", "Z", "PID"\n')
            f.write('ZONE T="P1", DATAPACKING=POINT, NODES={:d}, ELEMENTS={:d}, \
                    ZONETYPE=FETRIANGLE\n'.format(nvertex, nelement))
            # vertices
            for xyz in _xy:
                f.write("{:f} {:f} {:f} {:f}\n".format(xyz[0], xyz[1], xyz[2], 0.))
            # connectivity
            f.write("\n")
            for ijk in _edge:
                f.write("{:d} {:d} {:d}\n".format(int(ijk[0]), int(ijk[1]), int(ijk[2])))

"""
# Example
"""
if __name__ == "__main__":
    facet = Facets()
    # input file
    path_input = 'facets_vof.dat'
    facet.get_data_facets(path_input)

    # output file, in tecplot ASCII format
    path_output = 'facets_vof_tec.dat'
    facet.output_tecplot(path_output)