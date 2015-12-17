import random
import matplotlib.pyplot as plt
import numpy as np

class Person(object):
    def __init__(self, colour, location, racism=8):
        if colour in ["red", "blue"]:
            self.colour = colour
        elif colour == 1:
            self.colour = "red"
        elif colour == 2:
            self.colour = "blue"
        else:
            raise AttributeError('red or blue only')
        self.location = location
        self.racism = 10 - racism

    def dist(self, person):
        return np.sqrt((self.location[0]-person.location[0])**2+(self.location[1]-person.location[1])**2)

    def nearest(self, number, field):
        def get_key(item):
            return item[0]
        return np.array(sorted(np.array([[self.dist(person), person] for person in field.people]), key=get_key))[:number, 1]

    def happy(self, field):
        if len(filter(None, [self.colour != person.colour for person in self.nearest(10, field)])) > self.racism:
            return False
        else:
            return True

    def move(self, location):
        self.location = location

class Field(object):
    def __init__(self, n_persons):
        self.people = [Person(random.randint(1, 2), self.randloc()) for i in range(n_persons)]

    def randloc(self):
        return [10*np.random.standard_normal(), 10*np.random.standard_normal()]

    def show(self):
        x = [self.people[i].location[0] for i in range(len(self.people))]
        y = [self.people[i].location[1] for i in range(len(self.people))]
        colours = [self.people[i].colour for i in range(len(self.people))]
        sizes = [7**2 if person.happy(self) else 10**2 for person in self.people]
        plt.scatter(x, y, c=colours, s=sizes)
        plt.show()

    def move(self):
        for person in self.people:
            while not person.happy(self):
                person.move(self.randloc())

    def segregate(self, iters):
        for i in range(iters):
            field.show()
            field.move()
            if np.all([[person.happy(self) for person in self.people]]):
                field.show()
                break


field = Field(300)
field.segregate(10)
