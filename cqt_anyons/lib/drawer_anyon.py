import numpy as np


class DrawerAnyon:
    def __init__(self, initial_id, pos):
        self.__i = np.linspace(0, 1)
        self.__one = np.ones_like(self.__i)
        self.__x = [np.copy(self.__i)]
        self.__y = [pos * self.__one]

        self.__initial_id = initial_id

        self.__color = "black"
        self.__label = ""

    def __repr__(self):
        return f"{self.__label} (color: {self.__color}"

    def __str__(self):
        return f"{self.__label} (color: {self.__color}"

    @property
    def initial_id(self):
        return self.__initial_id

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, new):
        self.__label = new

    @property
    def color(self):
        return self.__color

    @color.setter
    def color(self, new):
        self.__color = new

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, new):
        self.__x.append(new)

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, new):
        self.__y.append(new)

    def get_last_x(self):
        return self.__x[-1][-1]

    def get_first_y(self):
        return self.__y[0][0]

    def get_last_y(self):
        return self.__y[-1][-1]

    def add_identity(self):
        self.__x.append(self.__i + self.get_last_x())
        self.__y.append(self.__one * self.get_last_y())
