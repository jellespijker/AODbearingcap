import matplotlib.pyplot as plt
from AOD.Material import *
from AOD.Model import *
from AOD.Bot import *
import warnings


def main():
    warnings.filterwarnings('ignore')
    model = Model()
    model.world.Layers['Soil'] = River_clay()
    [p_allow, p_load, depth, sink_depth, load] = model.solve_sinkdepth()
    torque_req = model.solve_torque(depth=sink_depth, load=load)
    print(torque_req.to('N*m'))
    print(sink_depth.to('mm'))
    print(load)
    plt.figure()
    plt.grid(True)
    plt.plot(load.item(0), -sink_depth.item(0), 'X')
    plt.plot(p_allow, -depth)
    plt.plot(p_load, -depth)
    plt.show()


if __name__ == '__main__':
    main()
