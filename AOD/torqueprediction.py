import matplotlib.pyplot as plt
from AOD.Material import *
from AOD.Model import *

from AOD.Bot import *


def main():
    # Define the fluid en solid layers
    air = Air()
    water = Water()
    silt = Silt()
    layers = {'Air': air, 'Fluid': water, 'Soil': silt}

    # define the world
    world = World(T=15.0 * ureg['degC'], layers=layers)

    # define the dredgebot
    dredgebot = Bot(no_of_screws=2)

    # define the model
    model = Model(world=world, bot=dredgebot)

    r = 1.
    D = np.arange(start=0, stop=r, step=0.1)
    theta = np.arccos((r - D) / r)
    B = 2 * r * np.sin(theta/2)

    plt.figure()
    plt.plot(D, B)
    plt.show()


if __name__ == '__main__':
    main()
