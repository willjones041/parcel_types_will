from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
class BaseParcel:
    def __init__(self, position):
        self.position = position
        

class PrecipitationParcel(BaseParcel):
    def __init__(self, position,  qe):
        super().__init__(position)
        self.qe = qe

class DynamicParcel(BaseParcel):
    def __init__(self, position,  ql):
        super().__init__(position)
        self.velocity = ql

class EulerianGrid:
    def __init__(self, grid_size):
                self.grid_size = grid_size
                self.grid = self.initialize_grid()

    def initialize_grid(self):
                return [[[None for _ in range(self.grid_size)] for _ in range(self.grid_size)] for _ in range(self.grid_size)]

    def plot_grid(self):
                

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                x, y, z = np.indices((self.grid_size, self.grid_size, self.grid_size))

                ax.scatter(x, y, z, c='b', marker='o')

                ax.set_xlabel('X axis')
                ax.set_ylabel('Y axis')
                ax.set_zlabel('Z axis')
                plt.show()
                
grid = EulerianGrid(grid_size=10)  # Specify the grid size
grid.plot_grid()